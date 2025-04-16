
import os,sys
sys.path.append(os.getcwd())
# from db.model_parser import PythiaParser
from db.Template.pythia_parser import PythiaParser
import re
import six
import scipy
import numpy as np

class HNLPythiaParser(PythiaParser):
    def __init__(self):
        super().__init__()
     
    def get_model_special_data(self):
        return self.make_interpolators()
        
    def make_interpolators(self, kind='linear'):
        """
        This function reads a file containing branching ratio histograms, and
        returns a dictionary of interpolators of the branching ratios, indexed by
        the decay string.
        """
        
        filepath = self.filegen.generate_filename("HNL", "Prod_brs","dat", name1 = "N1", name2 = "branchingratios" )
        histogram_data = self.parse_histograms(filepath)
        histograms = {}
        for (hist_string, (masses, br)) in six.iteritems(histogram_data):
            histograms[hist_string] = scipy.interpolate.interp1d(
                masses, br, kind=kind, bounds_error=False, fill_value=0, assume_sorted=True)
        return histograms

    def parse_histograms(self,filepath):
        """
        This function parses a file containing histograms of branching ratios.

        It places them in a dictionary indexed by the decay string (e.g. 'd_K0_e'),
        as a pair ([masses...], [branching ratios...]), where the mass is expressed
        in GeV.
        """
        with open(filepath, 'r') as f:
            lines = f.readlines()
        # Define regular expressions matching (sub-)headers and data lines
        th1f_exp      = re.compile(r'^TH1F\|.+')
        header_exp    = re.compile(r'^TH1F\|(.+?)\|B(?:R|F)/U2(.+?)\|.+? mass \(GeV\)\|?')
        subheader_exp = re.compile(r'^\s*?(\d+?),\s*(\d+?\.\d+?),\s*(\d+\.\d+)\s*$')
        data_exp      = re.compile(r'^\s*(\d+)\s*,\s*(\d+\.\d+)\s*$')
        # Locate beginning of each histogram
        header_line_idx = [i for i in range(len(lines)) if th1f_exp.match(lines[i]) is not None]
        # Iterate over histograms
        histograms = {}
        for offset in header_line_idx:
            # Parse header
            mh = header_exp.match(lines[offset])
            if mh is None or len(mh.groups()) != 2:
                raise ValueError("Malformed header encountered: {0}".format(lines[offset]))
            decay_code = mh.group(1)
            # Parse sub-header (min/max mass and number of points)
            ms = subheader_exp.match(lines[offset+1])
            if ms is None or len(ms.groups()) != 3:
                raise ValueError("Malformed sub-header encountered: {0}".format(lines[offset+1]))
            npoints  = int(ms.group(1))
            min_mass = float(ms.group(2))
            max_mass = float(ms.group(1))
            masses = np.linspace(min_mass, max_mass, npoints, endpoint=False)
            branching_ratios = np.zeros(npoints)
            # Now read the data lines (skipping the two header lines)
            for line in lines[offset+2:offset+npoints+1]:
                md = data_exp.match(line)
                if md is None or len(md.groups()) != 2:
                    raise ValueError("Malformed data row encountered: {0}".format(line))
                idx = int(md.group(1))
                br  = float(md.group(2))
                branching_ratios[idx] = br
            histograms[decay_code] = (masses, branching_ratios)
        return histograms
    
if __name__ == "__main__":
    test = HNLPythiaParser()
    
    print(test.get_model_special_data())
    