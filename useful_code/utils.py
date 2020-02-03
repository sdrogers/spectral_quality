import math
import numpy as np

def sqrt_normalise(peaks):
    temp = []
    total = 0.0
    for mz,intensity in peaks:
        temp.append((mz,math.sqrt(intensity)))
        total += intensity
    norm_facc = math.sqrt(total)
    normalised_peaks = []
    for mz,intensity in temp:
        normalised_peaks.append((mz,intensity/norm_facc))
    return normalised_peaks


class Spectrum(object):
    def __init__(self,peaks,file_name,scan_number,ms1,precursor_mz,parent_mz,rt = None,precursor_intensity = None,metadata = None):
        self.peaks = sorted(peaks,key = lambda x: x[0]) # ensure sorted by mz
        self.normalised_peaks = sqrt_normalise(self.peaks) # useful later
        self.n_peaks = len(self.peaks)
        self.max_ms2_intensity = max([intensity for mz,intensity in self.peaks])
        self.total_ms2_intensity = sum([intensity for mz,intensity in self.peaks])
        self.file_name = file_name
        self.scan_number = scan_number
        self.ms1 = ms1
        self.rt = rt
        self.precursor_mz = precursor_mz
        self.parent_mz = parent_mz
        self.precursor_intensity = precursor_intensity
        if metadata:
            self.metadata = metadata
        else:
            self.metadata = {}

    def get_annotation(self):
        if not 'annotation' in self.metadata:
            return None
        elif len(self.metadata['annotation']) == 0:
            return None
        else:
            anns = self.metadata['annotation']
            anns.sort(key = lambda x: x.score,reverse = True)
            return anns[0]



    annotation = property(get_annotation)


    def normalise_max_intensity(self,max_intensity = 1000.0):
        new_peaks = []
        for mz,intensity in self.peaks:
            new_peaks.append((mz,max_intensity*(intensity/self.max_ms2_intensity)))
        self.peaks = new_peaks
        self.n_peaks = len(self.peaks)
        if self.n_peaks > 0:
            self.normalised_peaks = sqrt_normalise(self.peaks)
            self.max_ms2_intensity = max([intensity for mz,intensity in self.peaks])
            self.total_ms2_intensity = sum([intensity for mz,intensity in self.peaks])
        else:
            self.normalised_peaks = []
            self.max_ms2_intensity = 0.0
            self.total_ms2_intensity = 0.0


    def randomise_intensities(self):
        from numpy.random import permutation
        intensities = [p[1] for p in self.peaks]
        permuted_intensities = permutation(intensities)
        new_peaks = []
        for i,(mz,intensity) in enumerate(self.peaks):
            new_peaks.append((mz,permuted_intensities[i]))
        self.peaks = new_peaks
        self.peaks.sort(key = lambda x: x[0])
        if self.n_peaks > 0:
            self.normalised_peaks = sqrt_normalise(self.peaks)
            self.max_ms2_intensity = max([intensity for mz,intensity in self.peaks])
            self.total_ms2_intensity = sum([intensity for mz,intensity in self.peaks])
        else:
            self.normalised_peaks = []
            self.max_ms2_intensity = 0.0
            self.total_ms2_intensity = 0.0        

    def flatten_peaks(self):
        max_intensity = max([p[1] for p in self.peaks])
        new_peaks = []
        for mz,intensity in self.peaks:
            new_peaks.append((mz,max_intensity))
        self.peaks = new_peaks
        if self.n_peaks > 0:
            self.normalised_peaks = sqrt_normalise(self.peaks)
            self.max_ms2_intensity = max([intensity for mz,intensity in self.peaks])
            self.total_ms2_intensity = sum([intensity for mz,intensity in self.peaks])
        else:
            self.normalised_peaks = []
            self.max_ms2_intensity = 0.0
            self.total_ms2_intensity = 0.0

    def remove_top_perc(self,perc):
        # remove the peaks corresponding to the top perc of intensity
        total_intensity = sum([p[1] for p in self.peaks])
        self.peaks.sort(key = lambda x: x[1])
        new_peaks = []
        total_found = 0.0
        for mz,intensity in self.peaks:
            total_found += intensity
            if total_found > (1-perc)*total_intensity:
                break
            else:
                new_peaks.append((mz,intensity))

        self.peaks = new_peaks
        self.peaks.sort(key = lambda x: x[0])

        if self.n_peaks > 0:
            self.normalised_peaks = sqrt_normalise(self.peaks)
            self.max_ms2_intensity = max([intensity for mz,intensity in self.peaks])
            self.total_ms2_intensity = sum([intensity for mz,intensity in self.peaks])
        else:
            self.normalised_peaks = []
            self.max_ms2_intensity = 0.0
            self.total_ms2_intensity = 0.0

    def remove_small_peaks(self,min_ms2_intensity = 10000):
        new_peaks = []
        for mz,intensity in self.peaks:
            if intensity >= min_ms2_intensity:
                new_peaks.append((mz,intensity))
        self.peaks = new_peaks
        self.n_peaks = len(self.peaks)
        if self.n_peaks > 0:
            self.normalised_peaks = sqrt_normalise(self.peaks)
            self.max_ms2_intensity = max([intensity for mz,intensity in self.peaks])
            self.total_ms2_intensity = sum([intensity for mz,intensity in self.peaks])
        else:
            self.normalised_peaks = []
            self.max_ms2_intensity = 0.0
            self.total_ms2_intensity = 0.0


    def remove_precursor_peak(self,tolerance = 17):
        new_peaks = []
        for mz,intensity in self.peaks:
            if abs(mz - self.precursor_mz) > tolerance:
                new_peaks.append((mz,intensity))
        self.peaks = new_peaks
        self.n_peaks = len(self.peaks)
        if self.n_peaks > 0:
            self.normalised_peaks = sqrt_normalise(self.peaks)
            self.max_ms2_intensity = max([intensity for mz,intensity in self.peaks])
            self.total_ms2_intensity = sum([intensity for mz,intensity in self.peaks])
        else:
            self.normalised_peaks = []
            self.max_ms2_intensity = 0.0
            self.total_ms2_intensity = 0.0
        

    def keep_top_k(self,k=6,mz_range=50):
        # only keep peaks that are in the top k in += mz_range
        start_pos = 0
        new_peaks = []
        for mz,intensity in self.peaks:
            while self.peaks[start_pos][0]< mz-mz_range:
                start_pos += 1
            end_pos = start_pos
            n_bigger = 0
            while end_pos < len(self.peaks) and self.peaks[end_pos][0] <= mz+mz_range:
                if self.peaks[end_pos][1] > intensity:
                    n_bigger += 1
                end_pos += 1
            if n_bigger < k:
                new_peaks.append((mz,intensity))

        self.peaks = new_peaks
        self.n_peaks = len(self.peaks)
        if self.n_peaks > 0:
            self.normalised_peaks = sqrt_normalise(self.peaks)
            self.max_ms2_intensity = max([intensity for mz,intensity in self.peaks])
            self.total_ms2_intensity = sum([intensity for mz,intensity in self.peaks])
        else:
            self.normalised_peaks = []
            self.max_ms2_intensity = 0.0
            self.total_ms2_intensity = 0.0
        

    def print_spectrum(self):
        print()
        print(self.file_name,self.scan_number)
        for i,(mz,intensity) in enumerate(self.peaks):
            print(i,mz,intensity,self.normalised_peaks[i][1])

    def plot(self,xlim = None,**kwargs):
        plot_spectrum(self.peaks,xlim=xlim,title = "{} {} (m/z= {})".format(self.file_name,self.scan_number,self.parent_mz),**kwargs)

    def __str__(self):
        return "Spectrum from scan {} in {} with {} peaks, max_ms2_intensity {}".format(self.scan_number,self.file_name,self.n_peaks,self.max_ms2_intensity)

    def __cmp__(self,other):
        if self.parent_mz >= other.parent_mz:
            return 1
        else:
            return -1

    def __lt__(self,other):
        if self.parent_mz <= other.parent_mz:
            return 1
        else:
            return 0


def make_spectrum(metadata,peaks):
    file_name = metadata['filename']
    scan_number = int(metadata['SCANS'])
    precursor_mz = float(metadata['PEPMASS'])
    try:
        rt = float(metadata['RTINSECONDS'])
    except:
        rt = None
    charge = metadata['CHARGE']
    ms1 = MS1(precursor_mz,rt,charge)
    return Spectrum(peaks,file_name,scan_number,ms1,precursor_mz,precursor_mz,rt=rt,metadata = metadata)


class MS1(object):
    def __init__(self,precursor_mz,rt,charge):
        self.precursor_mz = precursor_mz
        self.charge = charge
        self.rt = rt
        self.compute_parent_mz()
    def compute_parent_mz(self):
        try:
            self.parent_mz = self.precursor_mz*self.charge
            self.parent_mz -= PROTON_MASS*self.charge
        except:
            # if charge not available
            self.parent_mz = self.precursor_mz

def load_mgf(mgf_name,id_field = 'SCANS'):
    spectra  = {}
    with open(mgf_name,'r') as f:
        current_metadata = {'filename':mgf_name}
        current_peaks = []
        got_record = False
        for line in f:
            line = line.rstrip()
            if len(line) == 0:
                continue
            if line.startswith('BEGIN IONS'):
                if len(current_metadata) > 1:
                    if len(current_peaks) > 0:
                        spectrum = make_spectrum(current_metadata,current_peaks)
                        if id_field == 'SCANS':
                            id_val = int(current_metadata[id_field])
                        else:
                            id_val = current_metadata[id_field]
                        spectra[id_val] = spectrum
                        if len(spectra)%100 == 0:
                            print("Loaded {} spectra".format(len(spectra)))
                current_metadata = {'filename':mgf_name}
                current_peaks = []
            elif len(line.split('=')) > 1:
                # it is a metadata line
                tokens = line.split('=')
                current_metadata[tokens[0]] = tokens[1]
            elif not line.startswith('END IONS'):
                # it's a peak
                tokens = line.split()
                mz = float(tokens[0])
                intensity = float(tokens[1])
                current_peaks.append((mz,intensity))
    # save the last one
    if len(current_peaks) > 0:
        spectrum = make_spectrum(current_metadata,current_peaks)
        if id_field == 'SCANS':
            id_val = int(current_metadata[id_field])
        else:
            id_val = current_metadata[id_field]
        spectra[id_val] = spectrum
    return spectra


# Methods to compute the various stats from the proteomics paper

# number of peaks, sqrt transformed
def n_peaks_sqrt(spectrum):
    return math.sqrt(len(spectrum.peaks))

# arithmetic mean of the ion intensities, log transformed
def a_mean(spectrum):
    intensities = np.array([v[1] for v in spectrum.peaks])
    return np.log(intensities.mean())

def a_std(spectrum):
    intensities = np.array([v[1] for v in spectrum.peaks])
    return np.log(intensities.std())

def smallest_mz_range(spectrum,intensity_perc):
    # smallest mz range containing intensity_perc of the total intensity

    # sort by mz
    spectrum.peaks.sort(key = lambda x: x[0])

    intensities = np.array([v[1] for v in spectrum.peaks])    
    tot_intensity = intensities.sum()
    target_intensity = tot_intensity * intensity_perc
    
    smallest_mz_range = 1e10
    best_positions = (-1,-1)
    for i,peak in enumerate(spectrum.peaks):
        # check: is there still intensity_perc remaining?
        tot_left = intensities[i:].sum()
        if tot_left < target_intensity:
            break # we cannot get to the target starting here
        
        # otherwise, grow the region
        explained_intensity = intensities[i]
        pos = i
        while explained_intensity < target_intensity:
            pos += 1
            explained_intensity += intensities[pos]
        mz_range = spectrum.peaks[pos][0] - peak[0]
        if mz_range < smallest_mz_range:
            smallest_mz_range = mz_range
            best_positions = (i,pos)

    return smallest_mz_range,best_positions

def gap_std(spectrum):
    # std deviation of mz shifts, log transformed
    spectrum.peaks.sort(key = lambda x: x[0])
    mz = np.array([z[0] for z in spectrum.peaks])
    shifts = mz[1:] - mz[:-1]
    return np.log(shifts.std())

def n_neighbours(spectrum,width=2):
    spectrum.peaks.sort(key = lambda x: x[0])
    mz_counts = np.array([(z[0],0) for z in spectrum.peaks])
    for pos in range(len(spectrum.peaks) - 1):
        end_pos = pos + 1
        while end_pos < len(spectrum.peaks) and mz_counts[end_pos][0] < mz_counts[pos][0] + width:
            mz_counts[end_pos][1] += 1
            end_pos += 1
    for pos in range(len(spectrum.peaks)-1,0,-1):
        end_pos = pos - 1
        while end_pos >= 0 and mz_counts[end_pos][0] > mz_counts[pos][0] - width:
            mz_counts[end_pos][1] += 1
            end_pos -= 1
    n_neighbours = np.array([v[1] for v in mz_counts])
    return n_neighbours.mean()
    


