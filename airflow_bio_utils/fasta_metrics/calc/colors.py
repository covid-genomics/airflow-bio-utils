from .virdis_color_palette import VIRDIS_DATA
from matplotlib.colors import LinearSegmentedColormap
import matplotlib as mpl

VIRDIS_COLOR_MAP = LinearSegmentedColormap.from_list('viridis', VIRDIS_DATA)

LEGEND_MPL_BGCOLOR_PARAM_NAME = ''
from distutils.version import LooseVersion
if LooseVersion(mpl.__version__) <= '1.5.9':
    LEGEND_MPL_BGCOLOR_PARAM_NAME = 'axis_bgcolor'
else:
    LEGEND_MPL_BGCOLOR_PARAM_NAME = 'facecolor'