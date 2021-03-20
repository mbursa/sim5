# Raytracing class for disk models.
#
# Module defines a class that makes raytracing from disk photosphere
# based on supplied disk and spectral model.
# 
# This file is a part of SIM5 library. 
# See README and LICENCE file for details.


import sys
import math
import numpy as np
from sim5lib import * 


class DiskRaytrace:
    r_max = 1e6


    def __init__(self, bh_mass, bh_spin, bh_dist, disk_model, spectral_model):
        """
        Initializes the raytracing routine.

        Args:
            bh_mass: black hole mass [M_sun]
            bh_spin: black hole spin [0..0.999]
            bh_dist: observer distance [kpc]
            disk_model: (must be instance of <Sim5_DiskModel> class or a descendant)
            spectral_model:  (must be instance of <Sim5_DiskSpectrum> class or a descendant)
        """

        if (bh_spin < 1e-4): bh_spin = 1e-4
        
        self.bh_mass = bh_mass
        self.bh_spin = bh_spin
        self.bh_dist = bh_dist
        self.disk = disk_model
        self.spectra = spectral_model
    #end def



    def spectrum(self, incl, energies, limbdk=1, flat=0, radres=-1, angres=-1, hardening=0):
        """
        Computes disk spectrum.

        Args:
            incl: observer inclination [deg]
            energies: array of energies (at the detector) [keV]
            limbdk: limb darkening switch (limbdk>0 = on)
            radres: radial resolution factor (-1=use default)
            angres: angular resolution factor (-1=use default)
        Returns:
            observed spectrum [erg/s/cm2/keV]
        """
        #import sim5lib as sim5
        
        if (incl < 1.0): incl = 1.0
        
        of = open('bhspec.domega', 'w')
        of.truncate()
        
        spectrum1 = np.zeros(len(energies))
        spectrum2 = np.zeros(len(energies))
        spectrum_bb_f = np.zeros(len(energies))
        spectrum_bb_0 = np.zeros(len(energies))

        incl = math.radians(incl)
        if (radres<=0.0): radres = 0.15
        if (angres<=0.0): angres = 90.0

        nphi = int(math.floor(angres/math.sqrt(math.cos(incl))))
        dphi = 2.*math.pi/float(nphi);
        
        dOmega_tot = 0.0
        dOmega_err = 0.0

        rx = r_bh(self.bh_spin)
        while (rx<self.r_max*1.1):
            sys.stderr.write("Raytracing r=%.2f\n" % (rx))

            drx = radres*(1.+rx/5.)
            dOmega = math.cos(incl)*(rx+drx/2.)*drx*dphi * ((self.bh_mass*grav_radius)/(self.bh_dist*parsec*1e3))**2
            dOmega_tot += dOmega

            for iphi in range(nphi):
                phi = (float)(iphi) * dphi;
                alpha = -rx*math.cos(phi);
                beta  = -rx*math.sin(phi)*math.cos(incl)

                r, m, gd, k = self.geodesic(incl, alpha, beta, flat or (self.disk.h(1e5)==0.0))
                if (not gd): 
                    dOmega_err += dOmega
                    continue

                R = r*math.sqrt(1.-m*m)
                T = self.disk.t_eff(R)

                if (T == 0.0): continue

                tetrad = self.__tetrad(r,m)
                g = self.__gfactor(k, tetrad)
                e = self.__emission_angle(k, tetrad) if (limbdk>0) else -1.0

                if (not (e>=0.0 or e<1.0)): 
                    sys.stderr.write("invalid e r=%e e=%e\n"%(r,e))
                    #continue

                if (g == 0.0): 
                    sys.stderr.write("invalid g r=%e g=%e e=%e\n"%(r,g,e))
                    #continue

                if (not(g > 0.0)): continue

                S = self.disk.sigma(R)
                Q = self.__vertical_gravity(R, tetrad)
                f = hardening if (hardening>0) else self.__spectral_hardening(T, self.disk.lumi)

                #s, sok = self.spectra.spectrum(T, S, Q, e, f, energies/g)
                #spectrum1 += s*(g**3)*dOmega
                #if (sok): spectrum2 += s*(g**3)*dOmega
                spectrum_bb_f += self.spectra.spectrum(T, e, f, energies/g)*pow(g,3)*dOmega
                spectrum_bb_0 += self.spectra.spectrum(T, e, 1.0, energies/g)*pow(g,3)*dOmega
            #end of for

            rx = rx + drx
        #end of while
        sys.stderr.write("Raytracing completed.\n")

        of.close()

        #return spectrum1, spectrum2, spectrum_bb, spectrum_bb
        return spectrum_bb_f, spectrum_bb_0
    #end of def



    def image(self, incl, rmax, N, limbdk=1):
        """
        Computes disk image.

        Args:
            incl: observer inclination [deg]
            rmax: image extend [rg]
            N: pixels in image (both x and y dimension)
            limbdk: limb darkening switch (limbdk>0 = on)
        Returns:
            a tuple with image arrays for observed flux, gfactor, emission angle
        """
        # convert inclination to radians
        incl = math.radians(max(1.0,incl))

        # create image arrays
        image_F = np.full((N,N), None, dtype=np.float)
        image_g = np.full((N,N), None, dtype=np.float)
        image_e = np.full((N,N), None, dtype=np.float)
        image_R = np.full((N,N), None, dtype=np.float)
        image_T = np.full((N,N), None, dtype=np.float)
        image_H = np.full((N,N), None, dtype=np.float)
        image_V = np.full((N,N), None, dtype=np.float)

        sys.stderr.write("Imaging\n")
        for y in range(N):
            for x in range(N):
                alpha = ((x+.5)/N-0.5)*2.0*rmax
                beta  = ((y+.5)/N-0.5)*2.0*rmax
                dOmega = (2.0*rmax/N)**2 / ((self.bh_mass*grav_radius)/(self.bh_dist*parsec*1e3))**2
                # for log scale:
                #alpha = math.exp(alpha) if (alpha>0) else -math.exp(-alpha)
                #beta = math.exp(beta) if (beta>0) else -math.exp(-beta)

                r, m, gd, k = self.geodesic(incl, alpha, beta, flat=(self.disk.h(1e5)==0.0))

                if (not gd): continue
                
                R = r*math.sqrt(1.-m*m)
                F = self.disk.flux(R)
                T = self.disk.t_eff(R)
                V = self.disk.vr(R)

                if (F == 0.0): continue

                tetrad = self.__tetrad(r,m)
                g = self.__gfactor(k, tetrad)
                e = self.__emission_angle(k, tetrad)
                l = (0.5+0.75*e) if (limbdk>0) else 1.0

                if (e<0.0 or e>1.0): 
                    sys.stderr.write("invalid e r=%e e=%e\n"%(r,e))
                    if (e<0.0): e=0.0001
                    if (e>1.0): e=0.9999
                    #continue

                if (not(g > 0.0)): 
                    sys.stderr.write("invalid g (%e)\n" % (g))
                    continue
                
                image_F[y,x] = F * g**4 * l * dOmega
                image_g[y,x] = g
                image_e[y,x] = math.degrees(math.acos(e))
                image_T[y,x] = T
                image_R[y,x] = R
                image_H[y,x] = r*m
                image_V[y,x] = V
            #end for(x)
        #end for(y)
        sys.stderr.write("Imaging completed\n")

        return {'flux':image_F, 'gfactor':image_g, 'mue':image_e, 'T':image_T, 'R':image_R, 'H':image_H, 'V':image_V}
    #end of def



    def geodesic(self, incl, alpha, beta, flat=False):
        """
        Finds a geodesic with impact parameters [alpha,beta] that connects an observer
        at inclication `incl` and the disk surface.

        Args:
            incl: observer inclination [deg]
            alpha: impact parameter (horizontal) [GM/c2]
            beta: impact parameter (vertical) [GM/c2]
            flat: force flat disk (returns intersection with equatorial plane)
        Returns:
            a tuple of (radius, cos_theta, <geodesic> object, <doubleArray(4) object> with 4-momentum)
            or (0.0, 0.0, None) in case of error
        """
        status = intp()
        gd = geodesic()
        k = doubleArray(4)        
        r = m = 0.0

        geodesic_init_inf(incl, self.bh_spin, alpha, beta, gd, status)

        if (status.value() != 0):
            sys.stderr.write("not gd init\n")
            return (0.0, 0.0, None, None)

        if (flat):
            P = geodesic_find_midplane_crossing(gd, 0)
            r = geodesic_position_rad(gd, P)
            m = 0.0
        else:
            P, r, m = self.__find_surface(gd)
            if (not P): return (0.0, 0.0, None, None)
        #end of if

        if ((status.value()!=0) or math.isnan(r)): return (0.0, 0.0, None, None)

        photon_momentum(self.bh_spin, r, m, gd.l, gd.q, gd.Rpc-P, 1.0, k)
        
        return (r, m, gd, k)
    #end of def



    def __find_surface(self, gd, iteration=0):
        if (iteration>3): 
            #sys.stderr.write("__find_surface: too many iterations\n")
            return (None, 0, 0)

        # determine geometrically where the disk can intersect the ray
        disk_theta = math.atan(self.disk.h(1e6)/1e6)
        #r0 = math.sqrt( (gd.alpha**2+gd.beta**2)/(math.cos(gd.i)-tan_a*math.sin(gd.i))**2 * (1.+tan_a**2) )
        r0 = max(200.0, 1.1*gd.rp, (0.5+iteration)*math.sqrt((gd.alpha**2+gd.beta**2))/math.cos(gd.i+disk_theta))
        # set initial position for raytracing:
        # start at position that is sufficiently far and above the disk
        accuracy = 1e-2
        rbh = r_bh(self.bh_spin)

        while True:
            P1  = geodesic_P_int(gd, r0, 0)
            r1 = geodesic_position_rad(gd, P1)
            m1 = geodesic_position_pol(gd, P1)
            R1 = r1*math.sqrt(1.-m1*m1)
            H1 = r1*m1
            Hd = self.disk.h(R1)
            if ((Hd<H1) or (r0>5e6)): break;
            #sys.stderr.write("radius expand\n")
            r0 = 2.0*r0
        #end while

        # test if the initial point is above surface
        if (Hd >= H1): 
            #sys.stderr.write("__find_surface: abort bellow disk\n")
            return (None, 0, 0)

        #sys.stderr.write("r0=%e  hd=%e  z=%e\n"%(r0,Hd,H1))

        P = doublep(); P.assign(P1)
        r = doublep(); r.assign(r1)
        m = doublep(); m.assign(m1)
        status = intp()

        # follow trajectory with adaptive step until surface is reached
        step_factor = 1.0
        while True:
            # take step prop to sqrt(r)
            step = max(accuracy/2., min((H1-Hd)/2., 0.5*(math.sqrt(r.value())-0.99)*step_factor))
            geodesic_follow(gd, step, P, r, m, status)
            if (not status.value()): return (None, 0, 0)

            # update surface location
            R1 = r.value()*math.sqrt(1.-m.value()**2)
            H1 = r.value()*m.value()
            Hd = self.disk.h(R1)

            # surface hit?
            if (H1 <= Hd):
                if (step < accuracy):
                    geodesic_follow(gd, -step/2., P, r, m, status)
                    return (P.value(), r.value(), m.value())
                geodesic_follow(gd, -step, P, r, m, status)
                step_factor = step_factor/5.
                continue
            #/if

            # eq plane hit?
            if (H1 < 1e-4):
                P = geodesic_find_midplane_crossing(gd, 0)
                r = geodesic_position_rad(gd, P)
                m = geodesic_position_pol(gd, P)
                return (P, r, m)
            #/if

            # terminating conditions
            if (r.value() < 1.05*rbh): return (None, 0, 0)
            if (r.value() > 1.1*r0): return self.__find_surface(gd, iteration+1)
            if (m.value() < 0.0): return (None, 0, 0)
            if (step < accuracy/2.): break
        #/while

        # return error result
        #sys.stderr.write("__find_surface: other reason\n")
        return (None, 0, 0)
    #end of def



    def __tetrad(self, r, m):
        # return tetard attached to disk surface
        metric = sim5metric()
        tetrad = sim5tetrad()
        R = r*math.sqrt(1.-m*m)
        kerr_metric(self.bh_spin, r, m, metric)
        tetrad_surface(metric, Omega_from_ell(self.disk.l(R), metric),
            self.disk.vr(R), self.disk.dhdr(R) if m>0.0 else 0.0, tetrad)
        return tetrad
    #end of def



    def __gfactor(self, k, tetrad):
        # get energy shift (g-factor) with respect to infinity of photon that is INCOMING to the local position
        # note: since we raytrace backwards, we need to make the inversion of the momentum direction
        # note: 4-velocity is the first base vector of the tetrad and (1,0,0,0)*e[\mu][(a)] = U^\mu
        m = tetrad.metric
        U = doubleArray(4)
        on2bl(sim5vector((1,0,0,0)), U, tetrad)
        g = (double_array_getitem(k,0)*m.g00 + double_array_getitem(k,3)*m.g03) / dotprod(k, U, m)
        return g if (g>0.0) else 0.0
    #end of def



    def __gfactor_keplerian(self, r, l, q):
        if (r <= r_ms(self.bh_spin)): return 0.0
        OmegaK = 1./(self.bh_spin + math.pow(r,1.5))
        r2 = r**2
        a2 = self.bh_spin**2
        E = math.sqrt(1. - 2./r * (1.-self.bh_spin*OmegaK)**2 - (r2+a2)*(OmegaK**2));
        return E / (1. - OmegaK*l)
    #end of def



    def __emission_angle(self, k, tetrad):
        # cos(\mu_e) = (Z.N) = (P.N)/(P.U), where Z=U+P/(U.P) is unit space-like vector in direction of P,
        # N is surface normal vector and P is photon momentum vector
        m = tetrad.metric
        U = doubleArray(4)
        on2bl(sim5vector((1,0,0,0)), U, tetrad)
        N = doubleArray(4)
        on2bl(sim5vector((0,0,1,0)), N, tetrad)
        mue = dotprod(k, N, m)/dotprod(k, U, m)
        if math.isnan(mue): 
            sys.stderr.write("WRN: mue is nan (v1=%e, v2=%e)\n"%(dotprod(k, N, m), dotprod(k, U, m)))
        if (mue<0.0 and mue>-1e-2): mue = 1e-3  # due to imperfections in dH/dR derivative, mue can sometimes become slightly negative; we correct for that if it is not too big
        if (mue<0.0): sys.stderr.write("WRN: negative mue (%e, r=%e, m=%e)\n"%(mue, m.r, m.m))
        return mue
    #end of def



    def __emission_angle_keplerian(self, r, l, q):
        ar    = self.bh_spin*math.pow(r,-1.5)
        eol   = math.sqrt(1.-3.0/r+2.0*ar)/(1.+ar)
        kep   = 1./(math.pow(r,1.5)+self.bh_spin)
        mu_e  = abs(eol*math.sqrt(q)/(r*(1.-kep*l)))
        return mu_e
    #end of def



    def __vertical_gravity(self, R, tetrad):
        # get vertical gravity using formula from Zhu+2012
        # units: [s^-2]
        m = tetrad.metric
        U = doubleArray(4)
        on2bl(sim5vector((1,0,0,0)), U, tetrad)
        u_t = U[0]*m.g00 + U[3]*m.g03;
        u_f = U[0]*m.g03 + U[3]*m.g33;
        return self.bh_mass*solar_mass*grav_const/math.pow(R*self.bh_mass*grav_radius,3) * (u_f**2 + (self.bh_spin**2)*(u_t-1.))/R
    #end of def



    def __spectral_hardening(self, T, mdot):
        # returns disk hardening factor based on a fitting formula (You+2015; https://arxiv.org/abs/1506.03959).
        # T=effective temperature [K]; mdot=accretion rate [Mdot_Edd]
        t4 = T/1e4
        m4 = ((mdot+0.1)/0.2)**0.24
        if (t4 > 10.0):
            f = 1.6 * m4
        elif (t4 > 1.0):
            f = (t4/3)**0.3904 * m4
        else:
            f = m4
        return f
    #end of def



    def test_surface(self, incl):
        of = open('slimdisk_raytrace_surface.dat', 'w')
        of.truncate()

        incl = math.radians(incl)
        N = 50
        Rmax = 5.0

        for iy in range(N):
            for ix in range(N):
                alpha = (ix/N-0.5)*2*Rmax
                beta  = (iy/N-0.5)*2*Rmax
                print(ix, iy, alpha, beta)

                r, m, gd = self.geodesic(incl, alpha, beta, flat=False)
                if (not gd):
                    of.write("%d %d %e\n" % (iy, ix, 0))
                    continue
                else:
                    H = r*m
                    of.write("%d %d %e\n" % (iy, ix, H))
            #end of for
        #end of for

        of.close()
    #end of def


#end of class



