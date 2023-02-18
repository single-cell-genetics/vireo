from .version import __version__

from .utils import vcf_utils as vcf
from .utils import vireo_base as base
from .utils import vireo_model as model

from .utils.vcf_utils import load_VCF, match_SNPs
from .utils.io_utils import read_cellSNP, read_vartrix
from .utils.vireo_base import normalize, loglik_amplify, get_binom_coeff
from .utils.vireo_base import match, optimal_match

from .utils.vireo_wrap import vireo_wrap
from .utils.vireo_model import Vireo
from .utils.vireo_bulk import VireoBulk, LikRatio_test
from .utils.bmm_model import BinomMixtureVB

from . import plot

__all__ = [
    "__version__",
    "utils",
    "plot",
]