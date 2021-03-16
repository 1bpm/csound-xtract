#include "dtw.h"
#include <plugin.h>
#include "xtract/libxtract.h"
#include "xtract/xtract_stateful.h"
#include "xtract/xtract_scalar.h"
#include "xtract/xtract_helper.h"
#include <vector>
#include <stdexcept>
//#include "math.h"

#define MFCC_BANDS 13

MYFLT euclidean_distance(MYFLT* p1, MYFLT* p2, int bands) {
    MYFLT total = 0.0;
    for (int i = 0; i < bands; i++)
    {
        total = total + pow((p1[i] - p2[i]), 2);
    }
    return sqrt(total);
}

MYFLT manhattan_distance(MYFLT* p1, MYFLT* p2, int bands) {
    MYFLT sum = 0.0; 
    for (int i = 0; i < bands; i++) {
        for (int j = i + 1; j < bands; j++) {
            sum += (abs(p1[i] - p1[j]) + abs(p2[i] - p2[j])); 
        }
    }
    return sum; 
}




struct AnalysisProfile {
    bool centroid;
    bool mfccs;
    bool rms;
    bool zerocrossing;
	bool flatness;
	bool irregularity;
	bool power;
	bool sharpness;
	bool smoothness;
    int block_size;
    int buffer_size;
    
    int bands() {
        int num = 0;
        if (centroid) num ++;
        if (rms) num ++;
        if (zerocrossing) num ++;
		if (flatness) num ++;
		if (irregularity) num ++;
		if (power) num ++;
		if (sharpness) num ++;
		if (smoothness) num ++;
        if (mfccs) num += MFCC_BANDS;
        
        return num;
    }
    
    bool compatible(AnalysisProfile* t) {
        return (
            centroid == t->centroid
            && mfccs == t->mfccs
            && rms == t->rms
            && zerocrossing == t->zerocrossing
			&& flatness == t->flatness
			&& irregularity == t->irregularity
			&& power == t->power
			&& sharpness == t->sharpness
			&& smoothness == t->smoothness
        );
    }
};
	
struct AnalysisData {
    void* mutex;
    MYFLT* data;
    bool accumulated;
    int itemsPerBufferPeriod;
//    MYFLT timePerItem; 
    int position;
    bool ready;
    int sample_points;
    int data_size;
    AnalysisProfile* analysisProfile;

};

struct CorpusData {
    FUNC* input;
    AnalysisData* analysisData;
}; 


const char* badHandle = "cannot obtain data from handle";

template <typename T>
char* handleIdentifier() {
    if (std::is_same<T, AnalysisData>::value) { // TODO not working??
        return "::xtA%d";
    } else if (std::is_same<T, CorpusData>::value) {
        return "::xtC%d";
    } else if (std::is_same<T, AnalysisProfile>::value) {
        return "::xtP%d";
    }
    return "::xt%d";
}

/*
 * Obtain global object of typename from global variables by handle
 */
template <typename T>
T* getHandle(csnd::Csound* csound, MYFLT handle) {csnd::AuxMem<MYFLT> ax;
    char buffer[32];
    snprintf(buffer, 32, handleIdentifier<T>(), (int)handle);
    return (T*) csound->query_global_variable(buffer);  
}


/*
 * Create global object of typename in global variables, returning handle
 */
template <typename T>
MYFLT createHandle(csnd::Csound* csound, T** data) {
    char buffer[32];
    int handle = 0;
    snprintf(buffer, 32, handleIdentifier<T>(), handle);
    while ((*data = (T*) csound->query_global_variable(buffer)) != NULL) {
        snprintf(buffer, 32, handleIdentifier<T>(), ++handle);
    }
    csound->create_global_variable(buffer, sizeof(T));
    *data = (T*) csound->query_global_variable(buffer);
    
    return FL(handle);
}


class Analyser {
    csnd::Csound* csound;
    xtract_mel_filter mel_filters;
    MYFLT samplerate;
    MYFLT* windowed;
    MYFLT* spectrum;
    MYFLT* window;
	MYFLT* window_subframe;
    MYFLT* specargs;
	MYFLT* tempMfcc;
    int block_size;
    int buffer_size;
    int half_blocksize;
    AnalysisData* analysisOut;
    int outPos;
    bool firstIterationDone;
    bool accumulate;
    long sampling_length;
    int sample_point;
    int analysis_bands;
    
public:
    
    void init(
        csnd::Csound* csound, 
        AnalysisData* analysisOut, 
        bool accumulate=false,
        long sampling_length=0
    ) {
        
        this->csound = csound;
        this->sample_point = 0;
        this->sampling_length = sampling_length;
        this->accumulate = accumulate;
        this->analysisOut = analysisOut;

        analysis_bands = analysisOut->analysisProfile->bands();

        block_size = analysisOut->analysisProfile->block_size;
        buffer_size = analysisOut->analysisProfile->buffer_size;
        half_blocksize = block_size / 2;
        
        outPos = 0;
        analysisOut->ready = false;
        analysisOut->accumulated = accumulate;
        
        allocate_vectors();
        samplerate = csound->sr();
        windowed = (MYFLT*) csound->calloc(block_size * sizeof(MYFLT));
        spectrum = (MYFLT*) csound->calloc(block_size * sizeof(MYFLT));
        specargs = (MYFLT*) csound->calloc(4 * sizeof(MYFLT));
		tempMfcc = (MYFLT*) csound->calloc(MFCC_BANDS * sizeof(MYFLT));
        specargs[0] = samplerate / block_size;
        specargs[1] = XTRACT_MAGNITUDE_SPECTRUM; // XTRACT_LOG_POWER_SPECTRUM
        specargs[2] = 0; // DC component  
        specargs[3] = 0; // Normalisation  = 1
        
        window = xtract_init_window(block_size, XTRACT_HANN);
        window_subframe = xtract_init_window(block_size >> 1, XTRACT_HANN);
        xtract_init_wavelet_f0_state();
        
        if (analysisOut->analysisProfile->mfccs) {
            mel_filters.n_filters = MFCC_BANDS;
            mel_filters.filters = (MYFLT**) csound->malloc(MFCC_BANDS * sizeof(MYFLT*));

            for (int k = 0; k < MFCC_BANDS; k++) {
                mel_filters.filters[k] = (MYFLT*) csound->malloc(block_size * sizeof(MYFLT));
            }

            xtract_init_mfcc(
                block_size >> 1,    
                ((int)samplerate) >> 1, 
                XTRACT_EQUAL_GAIN, 
                20,
                20000, 
                mel_filters.n_filters, 
                mel_filters.filters
            );
        }
    }
    
    void allocate_vectors() {
        int items = 0;
        // TODO: some reasonable calculation here instead of loop
        for (int pos = 0; (pos + block_size) < buffer_size ; pos += half_blocksize) {
            items += 1;
        }
        analysisOut->itemsPerBufferPeriod = items;
        
        int sample_points = 1;
        if (accumulate) {
            sample_points = (int) (sampling_length / buffer_size);
        }
        analysisOut->sample_points = sample_points;
        analysisOut->data_size = sample_points * items * analysis_bands;
        analysisOut->data = (MYFLT*) csound->calloc(
            sizeof(MYFLT) * analysisOut->data_size
        );
        

    }
    
    void deinit() {
        //csound->free(windowed);
        //csound->free(spectrum);
        //csound->free(mfccargs);
        if (analysisOut->analysisProfile->mfccs) {
            for (int n = 0; n < MFCC_BANDS; ++n) {
                csound->free(mel_filters.filters[n]);
            }
            csound->free(mel_filters.filters);
        }
        xtract_free_window(window);
        xtract_free_window(window_subframe);
    }
     
    
    void analyse(MYFLT* buffer) {
        int buffer_period = 0;
        int data_pos;
        
        for (int pos = 0; (pos + block_size) < buffer_size ; pos += half_blocksize) {
            MYFLT *now = &buffer[pos];
            xtract_windowed(now, block_size, window, windowed);
            xtract_init_fft(block_size, XTRACT_SPECTRUM);
            xtract[XTRACT_SPECTRUM](windowed, block_size, &specargs[0], spectrum);
            xtract_free_fft();
            
            
            int data_pos;
            if (accumulate) {
                data_pos = analysisOut->itemsPerBufferPeriod * analysis_bands * sample_point + analysis_bands * buffer_period;
            } else {
                data_pos = analysisOut->itemsPerBufferPeriod * analysis_bands * 0 + analysis_bands * buffer_period;
            }
            
            if (analysisOut->analysisProfile->mfccs) {
                xtract_mfcc(spectrum, block_size >> 1, &mel_filters, tempMfcc);
				xtract_dct(tempMfcc, MFCC_BANDS, NULL, analysisOut->data+data_pos);
                data_pos += MFCC_BANDS;
            }

            if (analysisOut->analysisProfile->centroid) {
                xtract_spectral_centroid(spectrum, block_size, NULL, analysisOut->data+data_pos);
                data_pos ++;
            }
    
            if (analysisOut->analysisProfile->zerocrossing) {
                xtract_zcr(now, block_size, NULL, analysisOut->data+data_pos);
                data_pos ++;
            }

            if (analysisOut->analysisProfile->rms) {
                xtract_rms_amplitude(now, block_size, NULL, analysisOut->data+data_pos);
                data_pos ++;
            }
            
			if (analysisOut->analysisProfile->flatness) {
				xtract_flatness(spectrum, block_size >> 1, NULL, analysisOut->data+data_pos);
				data_pos ++;
			}
		
			if (analysisOut->analysisProfile->irregularity) {
				xtract_irregularity_j(spectrum, block_size >> 1, NULL, analysisOut->data+data_pos);
				data_pos ++;
			}

			if (analysisOut->analysisProfile->power) {
				xtract_power(spectrum, block_size >> 1, NULL, analysisOut->data+data_pos);
				data_pos ++;
			}

			if (analysisOut->analysisProfile->sharpness) {
				xtract_sharpness(spectrum, block_size >> 1, NULL, analysisOut->data+data_pos);
				data_pos ++;
			}

			if (analysisOut->analysisProfile->smoothness) {
				xtract_smoothness(spectrum, block_size >> 1, NULL, analysisOut->data+data_pos);
				data_pos ++;
			}

            
            buffer_period++;
        }
        if (accumulate) {
            sample_point++;
        }
    }
    
};


struct xtdump : csnd::Plugin<1, 1> {
    static constexpr char const *otypes = "k[]";
    static constexpr char const *itypes = "i";
    AnalysisData* input;
    int bands;
    
    int init() {
        if (!(input = getHandle<AnalysisData>(csound, inargs[0]))) {
            return csound->init_error("xtractor handle not valid");
        }
        bands = input->analysisProfile->bands();
        csnd::Vector<MYFLT> &out = outargs.vector_data<MYFLT>(0);
        out.init(csound, bands);
        return OK;
    }
    
    int kperf() {
        csnd::Vector<MYFLT> &out = outargs.vector_data<MYFLT>(0);       
        for (int i=0; i < bands; i++) {
            out[i] = *(input->data + i);
        }
        
        return OK;
    }
};

// kdone, kanalysis[] xtaccdump ixttractorhandle, koutputtrigger
struct xtaccdump : csnd::Plugin<2, 2> {
	static constexpr char const *otypes = "kk[]";
	static constexpr char const *itypes = "ik";
	AnalysisData* input;
	int bands;
	int accumulated;
	MYFLT* accumulation;

	int init() {
		if (!(input = getHandle<AnalysisData>(csound, inargs[0]))) {
			return csound->init_error("xtractor handle not valid");
		}
		bands = input->analysisProfile->bands();
		accumulated = 0;
		accumulation = (MYFLT*) csound->calloc(sizeof(MYFLT) * bands);
		outargs[0] = FL(0);
		csnd::Vector<MYFLT> &out = outargs.vector_data<MYFLT>(1);
		out.init(csound, bands);
		return OK;
	}

	void clearBuffer() {
		for (int i=0; i < bands; i++) {
			*(accumulation + i) = FL(0);
		}
	}

	int kperf() {
		csnd::Vector<MYFLT> &out = outargs.vector_data<MYFLT>(1);
		for (int i=0; i < bands; i++) {
			*(accumulation + i) += *(input->data + i);	
		}
		accumulated++;
		
		if (inargs[1] == FL(1)) {
			for (int i=0; i < bands; i++) {
				out[i] = *(accumulation + i) / accumulated;
			}
			outargs[0] = FL(1);
			accumulated = 0;
			clearBuffer();
		} else {
			outargs[0] = FL(0);
		}

		return OK;
	}
};


struct xtcorpusmatch : csnd::Plugin<1, 4> {
    static constexpr char const *otypes = "k";
    static constexpr char const *itypes = "iiko";
    CorpusData* corpus;
    AnalysisData* input;
    DTW::SimpleDTW evaluator;
	int bands;
    
    int init() {

        if (!(corpus = getHandle<CorpusData>(csound, inargs[0]))) {
            return csound->init_error("corpus handle not valid");
        }
        
        if (!corpus->analysisData->ready) {
            return csound->init_error("corpus analysis not ready");
        }
        
        if (!(input = getHandle<AnalysisData>(csound, inargs[1]))) {
            return csound->init_error("xtractor handle not valid");
        }

        if (!corpus->analysisData->analysisProfile->compatible(input->analysisProfile)) {
            return csound->init_error("xtractor profiles incompatible");
        }
	
		bands = corpus->analysisData->analysisProfile->bands();
        
        
        MYFLT (*distance_function)(MYFLT* p1, MYFLT* p2, int bands);
        switch ((int) inargs[3]) {
            case 0:
                distance_function = euclidean_distance;
                break;
            case 1:
                distance_function = manhattan_distance;
                break;
        }
        
        evaluator = DTW::SimpleDTW(
            csound,
            corpus->analysisData->itemsPerBufferPeriod, 
            input->itemsPerBufferPeriod,
            MFCC_BANDS,
            distance_function
        );
        
        return OK;
    }
    
    int nearest() {
        MYFLT best = INFINITY;
        MYFLT distance;
        int bestIndex = 0;
        
        
        for (int index = 0; index < corpus->analysisData->sample_points; index++) {
            distance = evaluator.EvaluateWarpingCost(
                input->data,
                input->itemsPerBufferPeriod, 
                corpus->analysisData->data+(corpus->analysisData->itemsPerBufferPeriod * bands * index + bands * 0),
                corpus->analysisData->itemsPerBufferPeriod
            );
            
            if (distance < best) {// && distance > 1) {
                best = distance;
                bestIndex = index;
            }
        } 
        return bestIndex * corpus->analysisData->analysisProfile->buffer_size;
    }
    
    int kperf() {
        if (inargs[2] == FL(1)) {
            outargs[0] = FL(nearest());
        } else {
            outargs[0] = 0;
        }
        return OK;
    }
};


struct xtprofile : csnd::Plugin<1, 11> {
    static constexpr char const *otypes = "i";
    static constexpr char const *itypes = "ooppppppppp";
    
    int init() {
        AnalysisProfile* profile;
        outargs[0] = createHandle<AnalysisProfile>(csound, &profile);
        profile->buffer_size = (inargs[0] != FL(0))? (int) inargs[0]: 4096;
        profile->block_size = (inargs[1] != FL(0))? (int) inargs[1]: 512;
        
        if (profile->block_size >= profile->buffer_size) {
            return csound->init_error("block size must be smaller than buffer size");
        }
        
        profile->mfccs = (bool) inargs[2];
        profile->centroid = (bool) inargs[3];
		profile->zerocrossing = (bool) inargs[4];
        profile->rms = (bool) inargs[5];
		profile->flatness = (bool) inargs[6];
		profile->irregularity = (bool) inargs[7];
		profile->power = (bool) inargs[8];
		profile->sharpness = (bool) inargs[9];
		profile->smoothness = (bool) inargs[10];
        return OK;
    }
};

/*
struct xtprofileprint : csnd::Plugin<1, 1> {
    static constexpr char const *otypes = "S";
    static constexpr char const *itypes = "i";
    int init() {
        AnalysisProfile* profile;
        if (!(profile = getHandle<AnalysisProfile>(csound, inargs[0]))) {
            return csound->init_error("profile handle invalid");
        }
        
        STRINGDAT* out = (STRINGDAT*) outargs(0);
        
        std::stringstream data;
        data << "\nBuffer size:\t" << profile->buffer_size
                << "\nBlock size:\t" << profile->block_size
                << "\nMFCCS:\t\t" << profile->mfccs
                << "\nCentroid:\t" << profile->centroid
                << "\nRMS:\t\t" << profile->rms << "\n";
        out->size = data.str().length();
        out->data = data.str().c_str();
        
        return OK;
    }
};
*/

// icorpushandle xtcorpus iprofilehandle, ifn 
struct xtcorpus : csnd::Plugin<1, 2> {
    static constexpr char const *otypes = "i";
    static constexpr char const *itypes = "ii";
    CorpusData* corpus;
    
    
    int init() {
        
        outargs[0] = createHandle<CorpusData>(csound, &corpus);
        
        AnalysisProfile* corpusProfile;
        if (!(corpusProfile = getHandle<AnalysisProfile>(csound, inargs[0]))) {
            return csound->init_error("profile handle invalid");
        }
       	CSOUND* csbase = (CSOUND *) csound->get_csound(); 
        if ((corpus->input = csbase->FTnp2Find(csbase, &inargs[1])) == NULL) {
            return csound->init_error("cannot get function table specified");
        }
        
        corpus->analysisData = (AnalysisData*) csound->malloc(sizeof(AnalysisData));
        corpus->analysisData->analysisProfile = corpusProfile;
        int buffer_size = corpusProfile->buffer_size;
        MYFLT* buffer;
        int buffer_position = 0;
        Analyser analyser;
        buffer = (MYFLT*) csound->malloc(sizeof(MYFLT) * buffer_size);
        analyser.init(csound, corpus->analysisData, true, corpus->input->flen);
        
        for (int i = 0; i < corpus->input->flen; i++) {
            buffer[buffer_position] = corpus->input->ftable[i];
            if (buffer_position == buffer_size -1) { 
                buffer_position = 0;
                analyser.analyse(buffer);
            } else {
                buffer_position++;
            }
        }

        //if (buffer_position != 0) {
        //    analyser.analyse(buffer, buffer_size-buffer_position);
        //}
        analyser.deinit();

        corpus->analysisData->ready = true;
        
        return OK;
    }
    
};



// ixthandle, kdone xtractor iprofile, ainput
struct xtractor : csnd::Plugin<2, 2> {
    static constexpr char const *otypes = "ik";
    static constexpr char const *itypes = "ia"; // trigger on input????
    Analyser analyser;
    AnalysisData* analysisData;
    AnalysisProfile* profile;
    int buffer_position;
    int ksmps;
    MYFLT* buffer;
    
    int init() {
        csound->plugin_deinit(this);
        buffer_position = 0;
        
        if (!(profile = getHandle<AnalysisProfile>(csound, inargs[0]))) {
            return csound->init_error("profile handle invalid");
        }
        
        buffer = (MYFLT*) csound->calloc(sizeof(MYFLT) * profile->buffer_size);
        ksmps = insdshead->ksmps;
        outargs[0] = createHandle<AnalysisData>(csound, &analysisData);

        analysisData->analysisProfile = profile;
        try {
            analyser.init(csound, analysisData, false);
        } catch (std::exception &exception) {
            return csound->init_error(exception.what());
        }
        return OK;
    }
    
    int deinit() {
        analyser.deinit();
        return OK;
    }
        
    int aperf() { 
        outargs[1] = FL(0);
        for (int i = 0; i < ksmps; i++) {
            buffer[buffer_position] = inargs(1)[i];
            if (buffer_position == profile->buffer_size -1) { 
                break;
            }
            buffer_position += 1;
        }
        
        if (buffer_position == profile->buffer_size - 1) {
            analyser.analyse(buffer);
            buffer_position = 0;
            outargs[1] = FL(1);
        }
        
        return OK;
    }
};

struct xtdistance : csnd::Plugin<1, 4> {
    static constexpr char const *otypes = "k";
    static constexpr char const *itypes = "iiko";
    AnalysisData* analysisData1;
    AnalysisData* analysisData2;
    DTW::SimpleDTW eval;
    
    int init() {
        if (!(analysisData1 = getHandle<AnalysisData>(csound, inargs[0]))) {
            return csound->init_error("xtractor handle 1 invalid");
        }
        
        if (!(analysisData2 = getHandle<AnalysisData>(csound, inargs[1]))) {
            return csound->init_error("xtractor handle 2 invalid");
        }
        
        if (!analysisData1->analysisProfile->compatible(analysisData2->analysisProfile)) {
            return csound->init_error("xtractor profiles incompatible");
        }
        
        MYFLT (*distance_function)(MYFLT* p1, MYFLT* p2, int bands);
        switch ((int) inargs[3]) {
            case 0:
                distance_function = euclidean_distance;
                break;
            case 1:
                distance_function = manhattan_distance;
                break;
        }
        
        eval = DTW::SimpleDTW(
            csound,
            analysisData1->itemsPerBufferPeriod,
            analysisData2->itemsPerBufferPeriod,
            MFCC_BANDS,
            distance_function
        );

        outargs[0] = INFINITY;
        return OK;
    }
    
    int kperf() {
        
        // && analysisData1->ready && analysisData1->ready
        if (inargs[2] == FL(1)) {
            outargs[0] = eval.EvaluateWarpingCost(
                analysisData1->data,
                analysisData1->itemsPerBufferPeriod,
                analysisData2->data,
                analysisData2->itemsPerBufferPeriod
            );
        } else {
            outargs[0] = INFINITY;
        }
         
        return OK;
    }
};







#include <modload.h>

void csnd::on_load(csnd::Csound *csound) {    
    csnd::plugin<xtprofile>(csound, "xtprofile", csnd::thread::i);
    //csnd::plugin<xtprofileprint>(csound, "xtprofileprint", csnd::thread::i);
    csnd::plugin<xtractor>(csound, "xtractor", csnd::thread::ia);
    csnd::plugin<xtdistance>(csound, "xtdistance", csnd::thread::ik);
    csnd::plugin<xtcorpus>(csound, "xtcorpus", csnd::thread::i);
    csnd::plugin<xtcorpusmatch>(csound, "xtcorpusmatch", csnd::thread::ik);
    csnd::plugin<xtdump>(csound, "xtdump", csnd::thread::ik);
	csnd::plugin<xtaccdump>(csound, "xtaccdump", csnd::thread::ik);
   
}
