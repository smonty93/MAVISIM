import numpy as np

class PSF:
    """PSF object class"""
    def __init__(self, data, metadata):
        self.data     = data
        self.dtype    = data.dtype
        self.xpos     = metadata["xpos"]
        self.ypos     = metadata["ypos"]
        self.pixsize  = metadata["pixsize"]
        self.Lambda   = metadata["Lambda"]
    # All other e2e specific psf methods can go here:

