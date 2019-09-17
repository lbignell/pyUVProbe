import numpy as np
import struct, time, os, re, string
import pandas as pd

#Data class will encapsulate the binary file info.
#Probably I'll add some analysis methods too.

class Data:
    #Read the binary file on initialisation.
    def __init__(self, fname, npoints):
        '''
        Read the UVProbe binary file fname.

        # Args:
        
        - fname: the name of the file (if not in current directory should be
          whole path).

        - npoints: Number of points in the UV-VIS spectrum. eg. for a spectrum
          taken between 200 and 700 nm at 0.5 nm intervals, there will be 901
          points. Must be integer.

        '''
        inputsOK = True
        self.binfile = fname
        self.npoints = npoints
        try:
            fh = open(fname, 'rb')
        except FileNotFoundError:
            print('File not found! Bailing...')
            inputsOK = False
       
        if not (type(npoints) is int):
            print('npoints must be integer, but type is {}! Bailing...'.format(type(npoints)))
            inputsOK = False

        if inputsOK:
            self._parsefile(fh)

    def _parsefile(self,fh):
        #The bits to read have been determined empirically, if something is
        #wrong after a version change, or a change of settings, these lines are
        #likely suspects.
        if self.npoints == 901:
            #Values for 250 -- 700 nm measurements:
            junk = fh.read(0x26c0)
            hdrblock = fh.read(0x2930 - 0x26c0)# Header info
            junk = fh.read(0x2a00 - 0x2930)
            firstblock = fh.read(0x4628 - 0x2a00) #contains the absorbance
            junk = fh.read(0x4800 - 0x4628)
            secondblock = fh.read(0x6428 - 0x4800)#contains the wavelengths
            fh.close()

        #Values for 300 -- 700 nm measurements:
        else:
            junk = fh.read(0x2680)
            hdrblock = fh.read(0x2950 - 0x2680)# Header info
            junk = fh.read(0x2a00 - 0x2950)
            firstblock = fh.read(self.npoints*8) #contains the absorbance
            junk = fh.read(30*8)
            secondblock = fh.read(self.npoints*8)#contains the wavelengths
            fh.close()
            #Getting data is easy, fortunately.
        
        try:
            self.abs = struct.unpack('<' + 'd'*self.npoints, firstblock)
            self.wl  = struct.unpack('<' + 'd'*self.npoints, secondblock)
        except:
            print('Error: len(firstblock) = {}, len(secondblock) = {}'.format(
                    len(firstblock), len(secondblock)))
            raise

        #Parsing the header is tricky. This may not work in the most general
        #case.
        self.txthdr = hdrblock.decode('cp850')#1252')
        #Hdr info will be stored in a dictionary
        self.hdrdir = {'Attachment Properties': {'Attachment': None}, 
                       'Instrument Properties': {'Instrument Type': None, 'Measuring Mode': None, 
                                                 'Slit Width': None, 'Accumulation time': None, 
                                                 'Light Source Change Wavelength': None,
                                                 'Detector Unit': None, 'S/R Exchange': None,
                                                 'Stair Correction': None}, 
                       'Sample Preparation Properties': {'Weight': None, 'Volume': None, 
                                                         'Dilution': None, 'Path Length': None,
                                                         'Additional Information': None}, 
                       'Measurement Properties': {'Wavelength Range (nm.)': None, 'Scan Speed': None,
                                                  'Sampling Interval': None, 
                                                  'Auto Sampling Interval': None, 'Scan Mode': None}}
        #Now fill the header info.
        for key, item in self.hdrdir.items():
            #Grab the chunk of string corresponding to the outer level header info
            startidx = self.txthdr.find('['+key+']') + len('[' + key + ']')
            tmphdr2 = self.txthdr[startidx:]
            endidx = tmphdr2.find('[')
            tmphdr3 = tmphdr2[:endidx]
            datainds_start = []
            datainds_end = []
            for innerkey, inneritem in item.items():
                datainds_start += [tmphdr3.find(innerkey)]
                datainds_end += [tmphdr3.find(innerkey) + len(innerkey)]
            datainds_start = np.array(datainds_start)
            datainds_end = np.array(datainds_end)
            for i, (innerkey, inneritem) in enumerate(item.items()):
            #Need to read between the end of this property and the start of the next
                diffs = datainds_start - datainds_end[i]
                diffs = diffs[diffs>0]
                if len(diffs) == 0:
                    tmphdr4 = tmphdr3[datainds_end[i]:]
                else:
                    tmphdr4 = tmphdr3[datainds_end[i]:(datainds_end[i] +  min(diffs))]
                item[innerkey] = ''.join([j for j in tmphdr4 if j in string.printable 
                                        and not j in string.whitespace])[1:]


