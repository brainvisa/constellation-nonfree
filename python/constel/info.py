version_major = 5
version_minor = 1
version_micro = 3
version_extra = ''

# Format expected by setup.py and doc/source/conf.py: string of form "X.Y.Z"
__version__ = "%s.%s.%s%s" % (version_major,
                              version_minor,
                              version_micro,
                              version_extra)

description = 'Constellation, Connectivity-based parcellation of the cortex'

long_description = """
=============
CONSTELLATION
=============

Connectivity-based parcellation of the cortex.

"""

# versions for dependencies
SPHINX_MIN_VERSION = '1.0'

# Main setup parameters
NAME = 'constellation-nonfree'
PROJECT = 'constellation'
ORGANISATION = "CEA"
MAINTAINER = "CEA"
MAINTAINER_EMAIL = ""
DESCRIPTION = description
LONG_DESCRIPTION = long_description
URL = ""
DOWNLOAD_URL = ""
LICENSE = ""
#CLASSIFIERS = CLASSIFIERS
AUTHOR = "CEA"
AUTHOR_EMAIL = ""
PLATFORMS = "Linux, Windows, Mac"
ISRELEASE = version_extra
VERSION = __version__
PROVIDES = ["constellation-gpl"]
REQUIRES = ["aims-free"]
EXTRAS_REQUIRE = {
    "doc": ["sphinx>=%s" % SPHINX_MIN_VERSION],
}
