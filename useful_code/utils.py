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