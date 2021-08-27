import vcf

from .segments_data_processor import SegmentsDataProcessor, SegmentsDataProcessor2
from .lst import lst
from .tai import tai
from .loh import loh


class HRD:
    """
    Class that preprocess input files and gets HRD scores for input sample
        
    Methods
    -------
    test_lst()
        Returns LST score/scores and DNA index of sample
        
    test_tai()
        Returns TAI score of sample
        
    test_loh()
        Returns LOH score of sample
        
    test_all()
        Returns LST, TAI, LOH scores, sum of these scores and DNA index of sample 
    """
    
    
    def __init__(self, seg_report_file, seg_report_file_with_header=True, seg_report_sample_name=None, vcf_file=None, vcf_sample_name=None):
        """
        Parameters
        ----------
        seg_report_file : str
            Path to segmental report
            
        seg_report_file_with_header: boolean
            Flag indicating if segmental report's file format has header
            
        seg_report_sample_name: str
            Name of sample in segmental report, if segmental report's file format doesn't have header
            
        vcf_file: str
            Path to VCF file. There must exist tabix file with same name and in same folder as provided VCF file.
        
        vcf_sample_name: str
            Name of sample in VCF file
        """
        
        if seg_report_file_with_header:
            self.sdp = SegmentsDataProcessor(seg_report_file)
        else:
            self.sdp = SegmentsDataProcessor2(seg_report_file, seg_report_sample_name)
            
        self.vcf_file = vcf_file
        self.vcf_sample_name = vcf_sample_name
        
        self.cnv_data = None
        self.ai_data = None
        self.loh_data = None
    
        
    def test_lst(self, LST_SMb=None):
        """
        Method that returns LST score and DNA index of input sample
        
        Parameters
        ----------
        LST_SMb=None: int, optional
            Value of parameter LST_SMb of LST method. If not provided, LST is count for value of parameter 3 - 11 Mb
            
        Returns
        -------
        Return value of lst function
        """
        
        if self.cnv_data is None:
            self.cnv_data = self.sdp.get_cnv_segments()
        vcf_reader = None
        if not self.vcf_file is None:
            vcf_reader = vcf.Reader(filename=self.vcf_file)
            
        return lst(self.cnv_data, vcf_reader, self.vcf_sample_name, LST_SMb)
    
    
    def test_tai(self):
        """
        Method that returns TAI score of input sample
        
        Returns
        -------
        Return value of tai function
        """
        
        if self.ai_data is None:
            self.ai_data = self.sdp.get_ai_segments()
        return tai(self.ai_data)
    
    
    def test_loh(self):
        """
        Method that returns LOH score of input sample
        
        Returns
        -------
        Return value of loh function
        """
        
        if self.loh_data is None:
            self.loh_data = self.sdp.get_loh_segments()
        return loh(self.loh_data)
    
    
    def test_all(self, LST_SMb=11):
        """
        Method that returns LST, TAI, LOH scores, sum of these scores and DNA index of input sample
        
        Parameters
        ----------
        LST_SMb=11: int, optional
            Value of parameter LST_SMb of LST method. If not provided, LST is count for value of parameter 3 - 11 Mb
        
        Returns
        -------
        scores: dict
            Dictionary with following structure:
            {
                LST: int
                    LST score of sample for parameter LST_SMb value
                TAI: int
                    TAI score of sample
                LOH: int
                    LOH score of sample
                HRD: int
                    Sum of LST, TAI and LOH scores of sample
                DNA index: float
                    DNA index of sample
            }
        """
        
        lst_score, dna_index = self.test_lst(LST_SMb)
        tai_score = self.test_tai()
        loh_score = self.test_loh()
        
        return {
            "LST": lst_score,
            "TAI": tai_score,
            "LOH": loh_score,
            "HRD": lst_score + tai_score + loh_score,
            "DNA index": dna_index
        }