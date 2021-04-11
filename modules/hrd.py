import vcf

from gap_data_processor import GapDataProcessor
from lengths_data_processor import process_lengths_data
from segments_data_processor import SegmentsDataProcessor, SegmentsDataProcessor2
from lst import lst
from tai import tai
from loh import loh


class HRD:
    
    def __init__(self, seg_report_file, seg_report_file_with_header=True, seg_report_sample_name=None, vcf_file=None, vcf_sample_name=None):
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
        if self.cnv_data is None:
            self.cnv_data = self.sdp.get_cnv_segments()
        vcf_reader = vcf.Reader(filename=self.vcf_file)
        return lst(self.cnv_data, vcf_reader, self.vcf_sample_name, LST_SMb)
    
    
    def test_tai(self):
        if self.ai_data is None:
            self.ai_data = self.sdp.get_ai_segments()
        return tai(self.ai_data)
    
    
    def test_loh(self):
        if self.loh_data is None:
            self.loh_data = self.sdp.get_loh_segments()
        return tai(self.loh_data)
    
    
    def test_all(self, LST_SMb=10):
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