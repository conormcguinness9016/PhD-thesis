#!/usr/bin/env python3
from SigProfilerExtractor import sigpro as sig
#sig.sigProfilerExtractor("vcf", "sigprofresultsSA", "germremoved_SA_SigProfileranalysis", "mm10", context_type="96", minimum_signatures=1, maximum_signatures=10, nmf_replicates=100)

sig.sigProfilerExtractor("vcf", "sigprofresultsSA_hets_sep", "germremoved_SA_SigProfileranalysis_Hets_sep", "mm10", context_type="96", minimum_signatures=1, maximum_signatures=10, nmf_replicates=100)
