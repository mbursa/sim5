import os
import sys
import math
import numpy as np
import importlib
import hashlib

class Sim5_FitsTable:

    def __init__(self, filename, ref_mass, ref_dist, params, energies, create_if_missing=True):
        """
        Creates in-memory FITS file to store spectral data.

        Args:
            ref_mass: reference BH mass [M_sun]
            ref_dist: reference BH distance [pc]
            params:   arrays of paramater grids [param1, param2, ...]
                         each param is a tuple with (grid_name, grid_values)
                         the latest param is changing fastest
            energies: energy grid [keV]
        """
        
        # import 'astropy.io.fits' module dynamicly to avoid package dependencies
        self.fits = importlib.import_module('astropy.io.fits')

        self.filename = filename
        print "fits:",filename
        self.params = params
        self.primary_hdu = None
        self.meta_hdu    = None
        self.spectra_hdu = None

        # calc params checksum
        m = hashlib.md5()
        m.update(str(ref_mass))
        m.update(str(ref_dist))
        for p in params: m.update(p[0]+str(p[1]))
        for e in energies: m.update(str(e))
        params_crc = m.hexdigest()
        sys.stderr.write("SlimULX_FitsWriter: CRC %s\n" % (params_crc))

        # calc total grid size
        self.total_grid_size = 1
        for i in range(len(params)): self.total_grid_size *= len(params[i][1])

        if (os.path.isfile(filename)):
            sys.stderr.write("SlimULX_FitsWriter: opening existing fits file %s\n" % (filename))
            try:
                f = self.fits.open(filename, memmap=False)  # open a FITS file
                self.primary_hdu = f[0]
                self.meta_hdu    = f[f.index_of('META')]
                self.spectra_hdu = f[f.index_of('SPECTRA')]
                if (params_crc != self.primary_hdu.header['CRC']): raise Exception('SlimULX_FitsWriter: cannot open, metadata differ')
            except:
                sys.stderr.write("SlimULX_FitsWriter: opening failed - invalid format\n")
                f.close()
                sys.exit()
            #end try
        #end if

        if (not self.meta_hdu):
            sys.stderr.write("SlimULX_FitsWriter: creating new fits file\n")

            # create primary HDU
            self.primary_hdu = self.fits.PrimaryHDU()
            self.primary_hdu.header['CRC'] = params_crc

            # create META table
            self.meta_hdu = self.fits.BinTableHDU.from_columns([
                self.fits.Column(name='NAME', format='16A'),
                self.fits.Column(name='N', format='1J'),
                self.fits.Column(name='GRID', format='1PE')
            ], None, len(params)+3)
            self.meta_hdu.name = 'META'
        
            # add reference values for mass and distance
            self.meta_hdu.data[0] = ('REF_MASS', 1, [ref_mass])
            self.meta_hdu.data[1] = ('REF_DIST', 1, [ref_dist])

            # add energy grid
            self.meta_hdu.data[2] = ('ENERGIES', len(energies), energies)

            # add grids
            for i in range(len(params)):
                param_name = params[i][0].upper()
                param_vals = params[i][1]
                self.meta_hdu.data[i+3] = (param_name, len(param_vals), param_vals)
            #end for
        
            # create spectral table
            self.spectra_hdu = self.fits.BinTableHDU.from_columns([
                #self.fits.Column(name='index', format='1J'),
                self.fits.Column(name='mdot', format='1E'),
                self.fits.Column(name='Iv_0', format=str(len(energies))+'E'),
                self.fits.Column(name='Iv_f', format=str(len(energies))+'E')
            ], None, self.total_grid_size)
            self.spectra_hdu.name = 'SPECTRA'

            self.save()
        #end if
        sys.stderr.write("SlimULX_FitsWriter: total grid size = %d\n" % (self.total_grid_size))
    #end def


    def generator(self):
        """
        Generator that iterates over the whole space of parameters.

        Because of the chosen format of the FITS file 
        the generator makes an outer loop over alpha_grid values and then an inner loop
        over spin+lumi+incl grid. It then yields a tuple containing an unique index for
        spectral data and values of alpha, spin, lumi and incl.

        Returns:
            Tuple with total index, individual paramster indexes and individual parameter values
        """
        index = 0
        while (index < self.total_grid_size):
            # skip row that have data already
            while (self.spectra_hdu.data[index][0] > 0.0): 
                index += 1
                if (index >= self.total_grid_size): break
            if (index >= self.total_grid_size): break

            # get indexes in each grid; the last grid is changing the fastest
            gindices = np.zeros(len(self.params), dtype=int)
            gvalues = np.zeros(len(self.params))
            N0 = self.total_grid_size
            for i in range(len(self.params)):
                N = len(self.params[i][1])
                N0 /= N
                gindices[i] = index//N0 % N
                gvalues[i] = self.params[i][1][gindices[i]]
    
            # return total index, array of grid indexes, and array of grid values
            #sys.stderr.write("> generetor: passing spectrum index %d (%s %s)\n" % (index, str(gindices), str(gvalues)))
            yield (index, gindices, gvalues)
            index += 1
        #end while
    #end def


    def write(self, index, mdot, Iv_0, Iv_f, flush=False):
        """
        Writes a spectrum at a position of `index`.
        """
        #self.spectra_hdu.data[index] = (index, spectrum)
        self.spectra_hdu.data[index] = [mdot, Iv_0, Iv_f]
        
        if flush: self.save()
    #end of def


    def save(self):
        """
        Writes the FITS table to a file.
        """
        sys.stderr.write("> saving FITS file")
        if os.path.isfile(self.filename): os.remove(self.filename)
        self.fits.HDUList([self.primary_hdu, self.meta_hdu, self.spectra_hdu]).writeto(self.filename)
        sys.stderr.write(" - done\n")
    #end of def

#end of class



