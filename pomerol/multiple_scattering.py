#!/usr/bin/python3
# -*-coding:utf-8 -*

import numpy as np
import matplotlib.pyplot as plt

from time import perf_counter #For the timer decorator

from atmosphere import *
from observation import *



DtoR = np.pi / 180.
RtoD = 1. / DtoR
RT = 6371. # km


class MultipleScattering:
    def __init__(self, in_dict, instrument_dir, atmosphere):

        self.atmosphere = atmosphere

        self.segment_max_length   = float(in_dict["MS_segment_max_length"]) #km
        self.segment_bins_length = float(in_dict["MS_segment_bins_length"]) #km

        self.min_altitude = float(in_dict["MS_min_altitude"]) #km
        self.max_altitude = float(in_dict["MS_max_altitude"]) #km

        self.max_events = int(in_dict['MS_max_events'])
        self.nb_rays = int(in_dict['MS_N_rays'])


        self.segment_distances = np.arange(0, self.segment_max_length, self.segment_bins_length)
        self.segment_Nbins_max = len(self.segment_distances)

        self.segment_mid_distances = self.segment_distances[:-1] + np.array([(self.segment_distances[i] - self.segment_distances[i-1]) / 2 for i in range(1, self.segment_Nbins_max)])

        # print(self.segment_mid_distances)

        self.initial_az = instrument_dir[0]
        self.initial_el = instrument_dir[1]
        self.initial_alt = instrument_dir[2]

        self.los = (np.sin(self.initial_el), np.cos(self.initial_el) * np.sin(self.initial_az), np.cos(self.initial_el) * np.cos(self.initial_az))

        self.hist = { "altitudes": [],
                      "elevations": [],
                      "lengths": [],
                      "sca_angles": [],
                      "nb_events": [],
                      "vectors": [],
                      "sca_planes": [],
                      "ray_cross_section": [],
                      "aer_cross_section": [],
                      "final_pos": [],
                      "total_length": []
        }

    # @timer
    def GetSegmentLength(self, altitude, elevation):
        """Compute the length of the path along one segment given its starting altitude and elevation.
        Compute atmospheric cross sections along the potential path which weigths the random probability to be scattered at one point.
        Return the segment length."""


        #compute the segment altitudes of each bin
        if elevation > 0:
            # segment_altitudes = np.arange(altitude, self.atmosphere.h_r_max, self.segment_bins_length * np.sin(elevation))
            segment_max_distance = self.segment_bins_length + (self.max_altitude - altitude) / np.sin(elevation) #+1 in case it passes the max altitude and escape to space!
        else:
            # segment_altitudes = np.arange(altitude, 0, self.segment_bins_length * np.sin(elevation))

            segment_max_distance = self.segment_bins_length + (self.min_altitude - altitude) / np.sin(elevation) #+1 in case it attains ground


        segment_max_bin = int(segment_max_distance / self.segment_bins_length)
        if segment_max_bin > 4999: #in case elevation ~0, the max number of bins of a list is 4999 ????!!!
            segment_max_bin = 4999
            segment_max_distance = 4999 * self.segment_bins_length
        segment_mid_distances = self.segment_mid_distances[:segment_max_bin]
        segment_altitudes = segment_mid_distances * np.sin(elevation)
        segment_altitudes += altitude

        # Get the cross section of one bin (ray+aer)
        cross_section = lambda alt: self.atmosphere.GetRSVolumeCS(alt) + self.atmosphere.GetAerosolCS(alt)

        # Get the cross section of each bin as defined above
        cross_sections = np.array([cross_section(a) for a in segment_altitudes])

        #For bin number n, compute the probability that scattering does not happen before on the segment (1-p(-1)) * ((1-p(-2))) * ... * ((1-p(-n)))
        # print(elevation*RtoD, altitude, segment_max_distance, segment_max_bin)
        cumul_anti_cs = np.zeros(segment_max_bin)
        cumul_anti_cs[0] = 1
        for i in range(1, segment_max_bin):
            cacs = (1 - cross_sections[i-1]) * cumul_anti_cs[i-1]
            cumul_anti_cs[i] = cacs
            # if cacs < 0.01:
            #     break
        cumul_anti_cs = np.array(cumul_anti_cs)

        proba = cross_sections * cumul_anti_cs
        # print(proba, sum(proba))
        proba = proba / sum(proba)

        # print(segment_altitudes)
        # print(cross_sections, sum(cross_sections))
        # print(cumul_anti_cs, sum(cumul_anti_cs))

        # f1, axs = plt.subplots(3, sharex = True, figsize=(16, 8))
        # # axs[0] = plt.subplot(111)
        # # axs[1] = plt.subplot(212)
        # # axs[2] = plt.subplot(313)
        #
        # axs[0].plot(segment_altitudes, cross_sections)
        # axs[1].plot(segment_altitudes, proba)
        # axs[2].plot(segment_altitudes, cumul_anti_cs)
        #
        # for a in axs:
        #     a.set_yscale('log')

        bin = np.random.choice(segment_max_bin, p = proba)

        length  = segment_mid_distances[bin]
        alt     = segment_altitudes[bin]
        cs_ray  = self.atmosphere.GetRSVolumeCS(alt)
        cs_aer  = self.atmosphere.GetAerosolCS(alt)

        return length, alt, cs_ray, cs_aer

    # @timer
    def GetScatteringAngle(self, cs_ray, cs_aer):
        phase_function = cs_ray * self.atmosphere.profiles["ray_Phase_Fct"] + cs_aer * self.atmosphere.profiles["aer_Phase_Fct"]
        phase_function /= (cs_ray + cs_aer)

        proba = phase_function / sum(phase_function)

        # f1, axs = plt.subplots(3, sharex = True, figsize=(16, 8))
        #
        # axs[0].plot(self.atmosphere.profiles["sca_angle"], proba)

        return np.random.choice(self.atmosphere.profiles["sca_angle"], p = proba)



    def GetScatteringPlane(self):
        ### Rotation matrix from instrument ref frame to UpEastNorth
        ca, sa = np.cos(self.initial_az), np.sin(self.initial_az)
        ce, se = np.cos(self.initial_el), np.sin(self.initial_el)
        Ria = np.array([	[se, 		0,		ce],
                            [	 ce * sa,	-ca,	-se * sa],
                            [	 ce * ca,	sa,		-se * ca]])

        plane_angle = np.random.random() * np.pi

        Pxi = 0
        Pyi = np.sin(plane_angle)
        Pzi = np.cos(plane_angle)

        Pu, Pe, Pn = np.dot(Ria, (Pxi, Pyi, Pzi))

        return np.cross(self.los, (Pu, Pe, Pn))

    # @timer
    def RotateAboutPlane(self, incident_vec, scattering_plane, angle):
        x, y, z = scattering_plane
        S = np.matrix([ [0, -z, y],
                        [z, 0, -x],
                        [-y, x, 0]])
        S2 = S @ S

        R = np.identity(3) + np.sin(angle) * S + (1-np.cos(angle)) * S2

        # print(R)
        # print(np.vstack(incident_vec))
        return R @ incident_vec

    # @timer
    def Propagate(self):
        N = self.max_events

        alt = self.initial_alt
        el  = self.initial_el
        scattering_plane = self.GetScatteringPlane()
        vec = np.vstack(self.los)

        alt_H = [alt]
        el_H  = [el]
        vec_H  = [vec]
        cs_ray_H  = []
        cs_aer_H  = []

        length_H = []
        sca_angle_H = []

        total_length = 0

        nb_events = 0
        for i in range(self.max_events):
            # print(f"MS Propagate: {i}, {alt}, {el*RtoD}, {scattering_plane}, {vec}")
            length, alt, cs_ray, cs_aer = self.GetSegmentLength(alt, el)
            nb_events += 1

            total_length += length

            if alt <= 0:
                # print(f"MS HIT GROUND after {i} scattering!!!")
                alt_H.append(alt)
                el_H.append(el)
                length_H.append(length)
                sca_angle_H.append(sca_angle)
                # N = nb_events+1
                break
            elif alt >= self.max_altitude:
                # print(f"MS ESCAPE to SPACE after {i} scattering!!!")
                alt_H.append(alt)
                el_H.append(el)
                length_H.append(length)
                sca_angle_H.append(sca_angle)
                # N = nb_events+1
                break

            sca_angle = self.GetScatteringAngle(cs_ray, cs_aer)
            # print(f"sca_angle: {sca_angle*RtoD}")

            length_H.append(length)
            sca_angle_H.append(sca_angle)

            vec = self.RotateAboutPlane(vec, scattering_plane, sca_angle)

            # print("vec", vec)
            el = np.arcsin(vec[0, 0])

            alt_H.append(alt)
            el_H.append(el)
            vec_H.append(vec)
            cs_ray_H.append(cs_ray)
            cs_aer_H.append(cs_aer)

        # print("DEBUG vec_H", len(vec_H), len(vec_H[0]))

        final_position = sum([v * l for v, l in zip(vec_H, length_H)])

        self.hist["altitudes"].append(alt_H)
        self.hist["elevations"].append(el_H)
        self.hist["lengths"].append(length_H)
        self.hist["sca_angles"].append(sca_angle_H)
        self.hist["vectors"].append(vec_H)
        self.hist["nb_events"].append(nb_events)
        self.hist["sca_planes"].append(scattering_plane)
        self.hist["ray_cross_section"].append(cs_ray_H)
        self.hist["aer_cross_section"].append(cs_aer_H)
        self.hist["final_pos"].append(final_position)
        self.hist["total_length"].append(total_length)

        # print("DEBUG MS.hist.vec_H", len(self.hist["vectors"]), len(self.hist["vectors"][0]), len(self.hist["vectors"][0][0]))

        # f_path, axs = plt.subplots(4, sharex=True)
        #
        # axs[0].plot(range(nb_events+1), alt_H)
        # axs[1].plot(range(nb_events+1), np.array(el_H) * RtoD)
        # axs[2].plot(range(nb_events), length_H)
        # axs[3].plot(range(nb_events), np.array(sca_angle_H)*RtoD)
        #
        # plt.show()

    @timer
    def PropagateAll(self):
        f_path, axs = plt.subplots(5, sharex=True)

        for i in range(self.nb_rays):
            self.Propagate()

            axs[0].plot(range(self.hist["nb_events"][i] + 1), self.hist["altitudes"][i], color="black", alpha=0.5)
            axs[1].plot(range(self.hist["nb_events"][i] + 1), np.array(self.hist["elevations"][i]) * RtoD, color="black", alpha=0.5)
            axs[2].plot(range(self.hist["nb_events"][i]), self.hist["lengths"][i], color="black", alpha=0.5)
            axs[3].plot(range(self.hist["nb_events"][i]), np.array(self.hist["sca_angles"][i])*RtoD, color="black", alpha=0.5)

        axs[4].hist(self.hist["nb_events"])


    def GetRayUnpolarisedFlux(self, ray_id, grd_map):

        u, e, n = self.hist["final_pos"][ray_id]
        print(u, e, n)

        az   = np.arctan2(e, n)
        dist = np.sqrt(e**2 + n**2)

        print(az, dist)

        F0 = grd_map.GetRadiantIntensity(az, dist)

        print("1", F0)
        if F0 == 0:
            return 0

        # F0 *= grd_map.GetArea(id)
        print("2", F0)
        print(self.hist["lengths"][ray_id])
        F0 /= np.prod([l**2 for l in self.hist["lengths"][ray_id] if l != 0])

        print("3", F0)

        return F0


    def GetTotalUnpolarisedFlux(self, grd_map):

        F = 0
        for ir in range(self.nb_rays):
            F += self.GetRayUnpolarisedFlux(ir, grd_map)

        print("F", F)
        return F



    def GetNextScatteringPos(self, init_pos, vec, length):
        new_pos = np.array(init_pos) + np.array(vec) * length
        return np.reshape(new_pos, 3)

    def GetRayPath(self, ray_id):
        pos = np.array([0, 0, 0])
        pos_H = [pos]

        # print(list(zip(self.hist["vectors"][ray_id], self.hist["lengths"][ray_id]))[0])

        for vec, length in zip(self.hist["vectors"][ray_id], self.hist["lengths"][ray_id]):
            # print(pos)

            next_pos = self.GetNextScatteringPos(pos, np.reshape(vec, 3), length)
            pos_H.append(next_pos)
            pos = next_pos

        return pos_H

    def Make3DPlot(self, grd_map = None):

        #Plot the path of all rays in 3D space
        f  = plt.figure()
        ax = plt.axes(projection ='3d')
        for id in range(self.nb_rays):
            pos_H = self.GetRayPath(id)
            u, e, n = zip(*pos_H)
            ax.plot(n, e, u, "*", alpha = 0.5)

        ### Plot the norm vector of all scattering planes in black and the line of sight in red
        f  = plt.figure()
        ax = plt.axes(projection ='3d')
        u, e, n = zip(*self.hist["sca_planes"])
        ax.plot(n, e, u, "k*", alpha = 0.5)
        ax.plot(self.los[2], self.los[1], self.los[0], "r*", alpha = 0.5)

    def MakeOriginPlot(self, grd_map):
        f, ax = grd_map.MakeInitialMapPlot()

        for iray in range(self.nb_rays):
            u, e, n = self.hist["final_pos"][iray]

            az   = np.arctan2(e, n)
            dist = np.sqrt(e**2 + n**2)

            ax.plot(az, dist, "r+")

    def MakeAltitudeHistogram(self):
        f, axs  = plt.subplots(1)

        altitudes_data = []
        for ray in self.hist["altitudes"]:
            altitudes_data.extend([a for a in ray if a != 0])

        w = np.ones_like(altitudes_data) / len(altitudes_data)
        hist, bins, _ = axs.hist(altitudes_data, weights=w)

        axs.set_title("MS Scattering events altitude histogram")

    def MakeScatteringHistograms(self, cs_ray, cs_aer):

        #Make histogram of N random scattering angle given the relative cross section of Rayleigh and Mie scaterring.
        N = 50000
        f, axs  = plt.subplots(1)

        sca_angles_data = [self.GetScatteringAngle(cs_ray, cs_aer)*RtoD for i in range(N)]

        w = np.ones_like(sca_angles_data) / len(sca_angles_data)
        hist, bins, _ = axs.hist(sca_angles_data, bins=self.atmosphere.profiles["sca_angle"]*RtoD, weights=w)

        phase_function = cs_ray * self.atmosphere.profiles["ray_Phase_Fct"] + cs_aer * self.atmosphere.profiles["aer_Phase_Fct"]
        phase_function /= cs_ray + cs_aer
        phase_function /= np.sum(phase_function)
        # phase_function /= 2

        print(np.sum(hist), np.sum(phase_function))

        # axs1 = axs.twinx()
        axs.plot(self.atmosphere.profiles["sca_angle"]*RtoD, phase_function, "r")

        axs.set_title("MS scattering angle histogram")
