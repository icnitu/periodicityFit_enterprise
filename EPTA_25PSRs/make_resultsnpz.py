import numpy as np
import pickle

#flagPSRs = ['J0751+1807', 'J0900-3144', 'J1012+5307', 'J1022+1001', 'J1600-3053', 'J1640+2224', 'J1713+0747', 'J1744-1134', 'J1909-3744', 'J1918-0642']
flagPSRs = ['J1713+0747']
for psr in flagPSRs:

    print(psr)

    with open('chains/'+psr+'/dynesty_results.pkl', 'rb') as f:
        data = pickle.load(f)

    np.savez('chains/'+psr+'/dynesty_samples.npz', data['samples'], data['logwt'], data['logz'])
