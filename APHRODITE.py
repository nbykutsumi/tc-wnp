import numpy as np
import sys, os, glob
from datetime import datetime
import myfunc.util as util
import calendar
import gzip
import struct

class v1101(object):
    def __init__(self, region="MA", res="050", dbbaseDir="/tank/utsumi/data/APHRO/APHRO_V1101"): 
        self.region=region
        self.res  = res
        self.dbbaseDir=dbbaseDir
        self.miss = -99.9
        if res=="050":
            dnynx = {"MA":[140,180],
                     "ME":[60,90],
                     "RU":[100,360],
                    }

            d1lat = {"MA": np.arange(-15+0.25, 55-0.25+0.001, 0.5),
                     "ME": np.arange( 15+0.25, 45-0.25+0.001, 0.5),
                     "RU": np.arange( 34+0.25, 84-0.25+0.001, 0.5),
                    }
            d1lon = {"MA": np.arange( 60+0.25, 150-0.25+0.001, 0.5),
                     "ME": np.arange( 20+0.25, 65 -0.25+0.001, 0.5),
                     "RU": np.arange( 15+0.25, 165-0.25+0.001, 0.5),
                    }
            d1latbnd = {"MA": np.arange(-15, 55+0.001, 0.5),
                        "ME": np.arange( 15, 45+0.001, 0.5),
                        "RU": np.arange( 34, 84+0.001, 0.5),
                       }
            d1lonbnd = {"MA": np.arange( 60, 150+0.001, 0.5),
                        "ME": np.arange( 20, 65 +0.001, 0.5),
                        "RU": np.arange( 15, 165+0.001, 0.5),
                       }


        elif res=="025":
            dnynx = {"MA":[280,360],
                     "ME":[120,180],
                     "RU":[200,720],
                    }
            d1lat = {"MA": np.arange(-15+0.125, 55-0.125+0.001, 0.25),
                     "ME": np.arange( 15+0.125, 45-0.125+0.001, 0.25),
                     "RU": np.arange( 34+0.125, 84-0.125+0.001, 0.25),
                    }
            d1lon = {"MA": np.arange( 60+0.125, 150-0.125+0.001, 0.25),
                     "ME": np.arange( 20+0.125, 65 -0.125+0.001, 0.25),
                     "RU": np.arange( 15+0.125, 165-0.125+0.001, 0.25),
                    }
            d1latbnd = {"MA": np.arange(-15, 55+0.001, 0.25),
                        "ME": np.arange( 15, 45+0.001, 0.25),
                        "RU": np.arange( 34, 84+0.001, 0.25),
                       }
            d1lonbnd = {"MA": np.arange( 60, 150+0.001, 0.25),
                        "ME": np.arange( 20, 65 +0.001, 0.25),
                        "RU": np.arange( 15, 165+0.001, 0.25),
                       }

        else:
            print("check res",res)
            sys.exit()

       
        self.ny, self.nx = dnynx[region]
        self.Lat = d1lat[region]
        self.Lon = d1lon[region]
        self.LatBnd = d1latbnd[region]
        self.LonBnd = d1lonbnd[region]
 
    def load_year(self, Year=None):
        res    = self.res
        region = self.region
        dbbaseDir = self.dbbaseDir 
        srcdir = dbbaseDir + "/APHRO_%s/%sdeg"%(region, res)
        ssearch= srcdir + "/APHRO_MA_050deg_V1101*.%04d.gz"%(Year)
        srcpath= glob.glob(ssearch)[0]
        ny = self.ny
        nx = self.nx
        if calendar.isleap(Year):
            nday = 366
        else:
            nday = 365

        with gzip.open(srcpath, "rb") as f:
            sdat = f.read()

        sfmt_all = "%df"%(ny*nx*2*nday)
        adat = np.array(struct.unpack(sfmt_all, sdat)).reshape(nday,2,ny,nx)

        return adat[:,0,:,:], adat[:,1,:,:]  # precip [mm/day], ratio_of_005box_with_stations [%]
                




