#!/usr/bin/python3
# -*-coding:utf-8 -*

import numpy as np


from mpi4py import MPI
mpi_comm = MPI.COMM_WORLD
mpi_rank = mpi_comm.Get_rank()
mpi_size = mpi_comm.Get_size()
mpi_name = mpi_comm.Get_name()

from time import perf_counter #For the timer decorator

from atmosphere import *
from observation import *
from simulation import *


DtoR = np.pi / 180.
RtoD = 1. / DtoR
RT = 6371. # km


class MultipleScattering:
    def __init__(self, in_dict, instrument_dir, atmosphere):

        self.atmosphere = atmosphere

        self.segment_max_length   = float(in_dict["MS_segment_max_length"]) #km
        self.segment_bins_length  = float(in_dict["MS_segment_bins_length"]) #km

        self.min_altitude = float(in_dict["MS_min_altitude"]) #km
        self.max_altitude = float(in_dict["MS_max_altitude"]) #km

        self.max_events = int(in_dict['MS_max_events'])
        self.nb_rays = int(in_dict['MS_N_rays'])

        self.nb_rays = int(np.ceil(self.nb_rays / mpi_size) * mpi_size)

        self.max_distance_from_instrument = float(in_dict["ground_emission_radius"]) #km

        self.segment_distances = np.arange(0, self.segment_max_length, self.segment_bins_length) #km

        self.segment_mid_distances = self.segment_distances[:-1] + self.segment_bins_length / 2. #np.array([(self.segment_distances[i+1] - self.segment_distances[i]) / 2 for i in range(0, self.segment_Nbins_max-1)])
        self.segment_Nbins_max = len(self.segment_mid_distances)

        # print(self.segment_mid_distances)

        self.initial_az = instrument_dir[0]
        self.initial_el = instrument_dir[1]
        self.initial_alt = instrument_dir[2]

        self.los = (np.sin(self.initial_el), np.cos(self.initial_el) * np.sin(self.initial_az), np.cos(self.initial_el) * np.cos(self.initial_az))

        self.hist = { "altitudes":[],
                      "elevations":[],
                      "lengths":[],
                      "sca_angles":[],
                      "nb_events":[],
                      "vectors":[],
                      "sca_planes":[],
                      "AoLPs":[],
                      "ray_cross_section":[],
                      "aer_cross_section":[],
                      "final_pos":[],
                      "total_length":[]
        }
        # self.hist = { "altitudes": [None] * self.nb_rays,
        #               "elevations": [None] * self.nb_rays,
        #               "lengths": [None] * self.nb_rays,
        #               "sca_angles": [None] * self.nb_rays,
        #               "nb_events": [None] * self.nb_rays,
        #               "vectors": [None] * self.nb_rays,
        #               "sca_planes": [None] * self.nb_rays,
        #               "AoLPs": [None] * self.nb_rays,
        #               "ray_cross_section": [None] * self.nb_rays,
        #               "aer_cross_section": [None] * self.nb_rays,
        #               "final_pos": [None] * self.nb_rays,
        #               "total_length": [None] * self.nb_rays
        # }

    # @timer
    def GetSegmentLength(self, altitude, elevation):
        """Compute the length of the path along one segment given its starting altitude and elevation.
        Compute atmospheric cross sections along the potential path which weigths the random probability to be scattered at one point.
        Return the segment length."""


        #Compute the distance to escape to space or touch ground
        if elevation >= 0:
            # segment_altitudes = np.arange(altitude, self.atmosphere.h_r_max, self.segment_bins_length * np.sin(elevation))
            segment_max_distance = (self.max_altitude - altitude) / np.sin(elevation)
            segment_end_altitude = self.max_altitude
        else:
            # segment_altitudes = np.arange(altitude, 0, self.segment_bins_length * np.sin(elevation))
            segment_max_distance = (self.min_altitude - altitude) / np.sin(elevation)
            segment_end_altitude = self.min_altitude


        segment_max_bin = int(segment_max_distance / self.segment_bins_length) #Number of bins in the segment
        # print("segment_max_bin", segment_max_bin)
        if segment_max_bin > self.segment_Nbins_max: #In case the segment is longer than very long (see input file). For example if elevation ~ 0
            segment_max_bin = self.segment_Nbins_max
            segment_max_distance = self.segment_max_length
            # print("segment_max_bin", segment_max_bin)

        segment_mid_distances = self.segment_mid_distances[:segment_max_bin]
        segment_altitudes = segment_mid_distances * np.sin(elevation)
        segment_altitudes += altitude

        # print("len", len(segment_mid_distances))

        # Get the volume cross section of one bin (ray+aer) in km-1
        cross_section = lambda alt: (self.atmosphere.GetRSVolumeCS(alt) + self.atmosphere.GetAerosolCS(alt))

        # Get the cross section of each bin as defined above km-1
        cross_sections = np.array([cross_section(a) for a in segment_altitudes])
        # print(cross_sections)

        # Compute the volume coefficient of the bin in km. I compute here only a proportionnal term because the proba is scaled to 1 in the end
        # volume_coeff = lambda l, dl: 3 * l**2 + (dl / 2)**2
        # volume_coeff = lambda l, dl: dl + dl**3 / (12 * l**2)
        # volumes = np.array([volume_coeff(l, self.segment_bins_length) for l in segment_mid_distances])

        cross_sections *= self.segment_bins_length

        # print(volumes, cross_sections)

        #Append one last bin to the lists that correspond to no scattering (either ground impact or escape to space), with cross section 1
        segment_max_bin += 1
        segment_mid_distances   = np.append(segment_mid_distances, segment_max_distance)
        segment_altitudes       = np.append(segment_altitudes, segment_end_altitude)
        cross_sections          = np.append(cross_sections, 1)

        # print(len(segment_mid_distances), len(cross_sections), segment_max_bin)

        #For bin number n, compute the probability that scattering does not happen before on the segment, i.e. probability that the ray reaches bin n. (1-p(-1)) * ((1-p(-2))) * ... * ((1-p(-n)))

        cumul_anti_cs = np.zeros_like(cross_sections) #Initialize the proba to zeros
        cumul_anti_cs[0] = 1 #First scattering bin along los: no scattering can happen before -> Probability 1 to reach first scattering bin

        for i in range(1, segment_max_bin):
            cacs = (1 - cross_sections[i-1]) * cumul_anti_cs[i-1]
            cumul_anti_cs[i] = cacs
            # if cacs < 0.01:
            #     break
        cumul_anti_cs = np.array(cumul_anti_cs)

        proba = cross_sections * cumul_anti_cs
        # print(proba, sum(proba))
        proba = proba / sum(proba)

        # if len(proba) < 10:
        #     print("alt", segment_altitudes)
            # print("proba", proba)

        # print(segment_altitudes)
        # print(cross_sections, sum(cross_sections))
        # print(cumul_anti_cs, sum(cumul_anti_cs))


        bin = np.random.choice(segment_max_bin, p = proba)
        end_ray = False
        if bin == segment_max_bin - 1:
            # print("Ray done !!")
            end_ray = True

        length  = segment_mid_distances[bin]
        alt     = segment_altitudes[bin]
        cs_ray  = self.atmosphere.GetRSVolumeCS(alt)
        cs_aer  = self.atmosphere.GetAerosolCS(alt)

        # f1, axs = plt.subplots(3, sharex = True, figsize=(16, 8))
        # # axs[0] = plt.subplot(111)
        # # axs[1] = plt.subplot(212)
        # # axs[2] = plt.subplot(313)
        #
        # axs[0].plot(segment_altitudes, cross_sections)
        # axs[1].plot(segment_altitudes, proba)
        # axs[2].plot(segment_altitudes, cumul_anti_cs)
        #
        #
        # for a in axs:
        #     a.set_yscale('log')

        return length, alt, cs_ray, cs_aer, end_ray

    # @timer
    def GetScatteringAngle(self, cs_ray, cs_aer):
        phase_function = cs_ray * self.atmosphere.profiles["ray_Phase_Fct"] + cs_aer * self.atmosphere.profiles["aer_Phase_Fct"]
        phase_function /= (cs_ray + cs_aer)

        proba = phase_function / sum(phase_function)

        # f1, axs = plt.subplots(3, sharex = True, figsize=(16, 8))
        #
        # axs[0].plot(self.atmosphere.profiles["sca_angle"], proba)

        random_angle = np.random.choice([-1, 1]) * np.random.choice(self.atmosphere.profiles["sca_angle"], p = proba)
        return random_angle


    def IsScattered(self, alt):
        cs_ray  = self.atmosphere.GetRSVolumeCS(alt)
        cs_aer  = self.atmosphere.GetAerosolCS(alt)

        cross_section = cs_ray + cs_aer
        cross_section *= self.segment_bins_length

        is_scattered = np.random.random() < cross_section

        if is_scattered:
            return True, cs_ray, cs_aer
        else:
            return False



    def GetScatteringPlane(self):
        ### Rotation matrix from instrument ref frame to UpEastNorth
        ca, sa = np.cos(self.initial_az), np.sin(self.initial_az)
        ce, se = np.cos(self.initial_el), np.sin(self.initial_el)
        Ria = np.array([	[se, 		0,		ce],
                            [	 ce * sa,	-ca,	-se * sa],
                            [	 ce * ca,	sa,		-se * ca]])

        ### Random plane angle
        plane_angle = np.random.random() * np.pi - np.pi/2

        ### Normal vector to plane in ref frame of instrument (==AoLP)
        Pxi = 0
        Pyi = np.sin(plane_angle)
        Pzi = np.cos(plane_angle)

        ### Normal vector to plane in UpEastNorth coord.
        Pu, Pe, Pn = np.dot(Ria, (Pxi, Pyi, Pzi))

        # return (0, 1, 0), plane_angle
        return (Pu, Pe, Pn), plane_angle

    # @timer
    def RotateAboutPlane(self, incident_vec, scattering_plane, angle):
        x, y, z = scattering_plane
        # S = np.matrix([ [0, -z, y],
        #                 [z, 0, -x],
        #                 [-y, x, 0]])
        # S2 = S @ S
        #
        # R = np.identity(3) + np.sin(angle) * S + (1-np.cos(angle)) * S2

        ca, sa = np.cos(angle), np.sin(angle)
        u = 1 - ca
        x2, y2, z2 = x**2, y**2, z**2

        R = np.matrix([ [ca + x2*u,         x*y*u - z*sa,       x*z*u + y*sa],
                        [x*y*u + z*sa,      ca + y2*u,          y*z*u - x*sa],
                        [x*z*u - y*sa,      y*z*u + x*sa,       ca + z2*u]])

        # print(R)
        # print(np.vstack(incident_vec))
        return R @ incident_vec



    # @timer
    def Propagate_old(self, id_ray = 0):
        N = self.max_events

        alt = self.initial_alt
        el  = self.initial_el
        scattering_plane, AoLP = self.GetScatteringPlane()

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
            length, alt, cs_ray, cs_aer, end_ray = self.GetSegmentLength(alt, el)
            nb_events += 1
            alt_H.append(alt)
            length_H.append(length)

            total_length += length

            if end_ray:
                break

            sca_angle = self.GetScatteringAngle(cs_ray, cs_aer)
            # print(f"sca_angle: {sca_angle*RtoD}")

            sca_angle_H.append(sca_angle)

            vec = self.RotateAboutPlane(vec, scattering_plane, sca_angle)

            # print("vec", vec)
            el = np.arcsin(vec[0, 0])

            el_H.append(el)
            vec_H.append(vec)
            cs_ray_H.append(cs_ray)
            cs_aer_H.append(cs_aer)

        # print("DEBUG vec_H", len(vec_H), len(vec_H[0]))

        #If the ray escape above a max altitude, do not record it!
        if alt_H[-1] == self.max_altitude:
            return False

        final_position = sum([v * l for v, l in zip(vec_H, length_H)])

        #If the ray ends up too far from the instrument (i.e. farther than the max radius of the ground map), do not record it!
        final_dist = np.sqrt(final_position[1]**2 + final_position[2]**2)
        if final_dist > self.max_distance_from_instrument:
            return False

        self.hist["altitudes"].append(alt_H)
        self.hist["elevations"].append(el_H)
        self.hist["lengths"].append(length_H)
        self.hist["sca_angles"].append(sca_angle_H)
        self.hist["vectors"].append(vec_H)
        self.hist["nb_events"].append(nb_events)
        self.hist["sca_planes"].append(scattering_plane)
        self.hist["AoLPs"].append(AoLP)
        self.hist["ray_cross_section"].append(cs_ray_H)
        self.hist["aer_cross_section"].append(cs_aer_H)
        self.hist["final_pos"].append(final_position)
        self.hist["total_length"].append(total_length)

        # self.hist["altitudes"][id_ray] = alt_H
        # self.hist["elevations"][id_ray] = el_H
        # self.hist["lengths"][id_ray] = length_H
        # self.hist["sca_angles"][id_ray] = sca_angle_H
        # self.hist["vectors"][id_ray] = vec_H
        # self.hist["nb_events"][id_ray] = nb_events
        # self.hist["sca_planes"][id_ray] = scattering_plane
        # self.hist["AoLPs"][id_ray] = AoLP
        # self.hist["ray_cross_section"][id_ray] = cs_ray_H
        # self.hist["aer_cross_section"][id_ray] = cs_aer_H
        # self.hist["final_pos"][id_ray] = final_position
        # self.hist["total_length"][id_ray] = total_length

        # print("DEBUG MS.hist.vec_H", len(self.hist["vectors"]), len(self.hist["vectors"][0]), len(self.hist["vectors"][0][0]))

        # f_path, axs = plt.subplots(4, sharex=True)
        #
        # axs[0].plot(range(nb_events+1), alt_H)
        # axs[1].plot(range(nb_events+1), np.array(el_H) * RtoD)
        # axs[2].plot(range(nb_events), length_H)
        # axs[3].plot(range(nb_events), np.array(sca_angle_H)*RtoD)
        #
        # plt.show()

        return True



    def Propagate(self, alt, vec, scattering_plane, id_ray = 0, hist = np.zeros((0, 6)), segment_length = 0):

        is_scattered = False

        # While the ray is not scattered == while the ray goes in a straight line == while we are on the same path segment
        while not is_scattered:

            # print(alt, vec, vec[0, 0], segment_length, self.segment_bins_length)

            alt += (vec[0, 0] * self.segment_bins_length) #Increment the altitude
            segment_length += self.segment_bins_length #Increment the segment length

            # print(alt, vec, vec[0, 0], segment_length, self.segment_bins_length)

            if alt < self.min_altitude:
                hist = np.append(hist, [[alt, None, None, None, None, segment_length]], axis = 0)
                return True, hist
            elif alt > self.max_altitude:
                return False, hist

            is_scattered = self.IsScattered(alt)

        _, cs_ray, cs_aer = is_scattered
        sca_angle = self.GetScatteringAngle(cs_ray, cs_aer)
        # print(vec, sca_angle*RtoD)
        vec = self.RotateAboutPlane(vec, scattering_plane, sca_angle)
        # print(vec)

        hist = np.append(hist, [[alt, vec, sca_angle, cs_ray, cs_aer, segment_length]], axis = 0)

        # print(hist)

        segment_length = 0

        return self.Propagate(alt, vec, scattering_plane, id_ray=id_ray, hist = hist, segment_length=segment_length)


    @timer
    def PropagateAll(self):
        # f_path, axs = plt.subplots(5, sharex=True)

        for id_ray in range(int(self.nb_rays / mpi_size)):
            # id_ray = i * mpi_size + mpi_rank
            # print(i, id_ray, self.nb_rays)

            alt = self.initial_alt
            scattering_plane, AoLP = self.GetScatteringPlane()

            vec = np.vstack(self.los)


            invalid_ray = True
            while invalid_ray:
                hist = np.array([[alt, vec, None, None, None, None]])
                # print("START OVER")
                hit_ground, hist = self.Propagate(alt, vec, scattering_plane, id_ray = id_ray, hist = hist)
                if not hit_ground:
                    continue

                vectors         = hist[:-1, 1]
                lengths         = hist[1:, 5]
                final_position  = sum(vectors * lengths)

                if np.linalg.norm(final_position) > self.max_distance_from_instrument:
                    continue

                invalid_ray = False


            nb_events = hist.shape[0] - 2 # Number of scattering events (not counting the start and end of the ray)
            # print("nb_events", nb_events)
            altitudes   = hist[:, 0]
            # print("altitudes", altitudes)
            # print("vectors", vectors)
            elevations  = np.arcsin([v[0, 0] for v in vectors])
            # print("elevations", elevations*RtoD)
            sca_angles  = hist[1:-1, 2]
            # print("sca_angles", sca_angles*RtoD)
            # print("lengths", lengths)

            total_length    = sum(lengths)
            # print("total_length", total_length)

            # print("final_position", final_position)

            ray_cross_section = hist[1:-1, 3]
            # print("ray_cross_section", ray_cross_section)
            aer_cross_section = hist[1:-1, 4]
            # print("aer_cross_section", aer_cross_section)

            self.hist["altitudes"].append(altitudes)
            self.hist["elevations"].append(elevations)
            self.hist["lengths"].append(lengths)
            self.hist["sca_angles"].append(sca_angles)
            self.hist["vectors"].append(vectors)
            self.hist["nb_events"].append(nb_events)
            self.hist["sca_planes"].append(scattering_plane)
            self.hist["AoLPs"].append(AoLP)
            self.hist["ray_cross_section"].append(ray_cross_section)
            self.hist["aer_cross_section"].append(aer_cross_section)
            self.hist["final_pos"].append(final_position)
            self.hist["total_length"].append(total_length)

            # print(self.hist)



        #     axs[0].plot(range(self.hist["nb_events"][i] + 1), self.hist["altitudes"][i], color="black", alpha=0.5)
        #     axs[1].plot(range(self.hist["nb_events"][i]), np.array(self.hist["elevations"][i]) * RtoD, color="black", alpha=0.5)
        #     axs[2].plot(range(self.hist["nb_events"][i]), self.hist["lengths"][i], color="black", alpha=0.5)
        #     axs[3].plot(range(self.hist["nb_events"][i]-1), np.array(self.hist["sca_angles"][i])*RtoD, color="black", alpha=0.5)
        #
        # axs[4].hist(self.hist["nb_events"])

        mpi_comm.Barrier()
        for k in self.hist.keys():
            # if mpi_rank == 0: print(k, len(self.hist[k]))
            receiver = mpi_comm.gather(self.hist[k][:], root=0)
            # if mpi_rank == 0: print(len(self.hist[k]), self.hist[k])
            if mpi_rank == 0:
                for i in range(1, mpi_size):
                    self.hist[k].extend(receiver[i])
            # if mpi_rank == 0: print(len(self.hist[k]), self.hist[k])



    def GetTotalContribution(self):

        # print(self.Flux_list, self.Stocks_I_list, self.Stocks_Q_list)
        V = 0
        Vc = 0
        Vs = 0

        i_l = []
        d_l = []
        a_l = []

        for i_ray in range(self.nb_rays):
            I    = self.Stocks_I_list[i_ray]
            DoLP = self.Stocks_Q_list[i_ray]
            AoLP = self.hist["AoLPs"][i_ray]

            if DoLP < 0:
                DoLP *= -1
                AoLP += np.pi/2

            V += I
            Vc += I * DoLP * np.cos(2 * AoLP) / 2
            Vs += I * DoLP * np.sin(2 * AoLP) / 2

            DoLP = 100 * 2 * np.sqrt(Vc**2 + Vs**2) / V
            AoLP = np.arctan2(Vs, Vc) / 2.
            # print(DoLP, AoLP*RtoD)

            i_l.append(V)
            d_l.append(DoLP)
            a_l.append(AoLP*RtoD)

            # print(DoLP, AoLP*RtoD)

        # DoLP = 100 * np.sqrt(Q**2 + U**2) / I
        # AoLP = np.arctan2(U, Q) / 2.
        DoLP = 100 * 2 * np.sqrt(Vc**2 + Vs**2) / V
        AoLP = np.arctan2(Vs, Vc) / 2.

        f, axs = plt.subplots(3, sharex=True)
        axs[0].plot(range(self.nb_rays), i_l, "k-")
        axs[1].plot(range(self.nb_rays), d_l, "r-")
        axs[2].plot(range(self.nb_rays), a_l, "b-")

        # print(V, Vc, Vs, DoLP, AoLP*RtoD)
        return V, Vc, Vs, DoLP, AoLP

        # return S_I_total, S_Q_total

    def SetStocksParameters(self):
        """Compute the stocks paramters of each individual rays and stock them in a list"""
        self.Stocks_I_list = np.zeros(self.nb_rays)
        self.Stocks_Q_list = np.zeros(self.nb_rays)

        for i_ray in range(self.nb_rays):
            self.Stocks_I_list[i_ray], self.Stocks_Q_list[i_ray] = self.GetRayStocksParameters(i_ray)

            self.Stocks_I_list[i_ray] *= self.Flux_list[i_ray]

        # print("self.Stocks_I_list", self.Stocks_I_list)
        # print("self.Stocks_Q_list", self.Stocks_Q_list)

    def GetRayStocksParameters(self, ray_id):
        """Forward scattering of a single ray. Th ray is supposed emited unpolarized (1, 0)."""
        S = np.vstack((1, 0))

        for i_sca in range(self.hist["nb_events"][ray_id], 0, -1):
            theta = self.hist["sca_angles"][ray_id][i_sca - 1]
            M = self.GetTotalScatteringMatrix(theta, self.hist["ray_cross_section"][ray_id][i_sca - 1], self.hist["aer_cross_section"][ray_id][i_sca - 1])

            S = M @ S
        # print("ray stocks", S)
        return S


    def GetRayleighScatteringMatrix(self, a):

        ca = np.cos(a)
        M11 = ca**2 + 1
        M12 = ca**2 - 1

        M = 3/4 * np.matrix([[M11, M12], [M12, M11]])

        return M

    def GetMieScatteringMatrix(self, a):
        """NOT SURE ABOUT THAT!!!"""
        i, d = self.atmosphere.GetAerosolPhaseFunction(a)

        M = np.matrix([[i, 0], [0, d]])

        # return M
        return 0

    def GetTotalScatteringMatrix(self, a, ray_cs, aer_cs):
        """Average of the rayleigh and Mie scattering scattering matrices."""
        M_rs = self.GetRayleighScatteringMatrix(a)
        M_mie = self.GetMieScatteringMatrix(a)

        M = ray_cs * M_rs + aer_cs * M_mie
        return M / (ray_cs + aer_cs)


    def GetTotalUnpolarisedFlux(self, grd_map):

        F = np.sum(self.Flux_list)

        # print("F", F)
        return F

    def SetRaysFluxList(self, grd_map):
        self.Flux_list = np.zeros(self.nb_rays)
        for ir in range(self.nb_rays):
            self.Flux_list[ir] = self.GetRayUnpolarisedFlux(ir, grd_map)


    def GetRayUnpolarisedFlux(self, ray_id, grd_map):

        # print(self.hist["final_pos"][ray_id])
        u, e, n = self.hist["final_pos"][ray_id]
        # print(u, e, n)

        az   = np.arctan2(e, n)
        dist = np.sqrt(e**2 + n**2)

        # print(az*RtoD, dist)

        F0 = grd_map.GetRadiantIntensity(az, dist)

        # print("1", F0)
        if F0 == 0:
            return 0

        # F0 *= grd_map.GetArea(id)
        # print("2", F0)
        # print(self.hist["lengths"][ray_id])
        F0 /= np.prod([(l*1000)**2 for l in self.hist["lengths"][ray_id] if l != 0])

        # print("3", F0)

        return F0




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
            ax.plot([self.los[2], self.los[2]* 50], [self.los[1], self.los[1]* 50], [self.los[0], self.los[0]* 50], "r-", alpha = 0.5)

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
            altitudes_data.extend([a for a in ray if a > 0])

        w = np.ones_like(altitudes_data) / len(altitudes_data)
        hist, bins, _ = axs.hist(altitudes_data, bins = 100, weights=w)

        axs.set_title("MS Scattering events altitude histogram")

    def MakeScatteringHistograms(self, cs_ray=0, cs_aer=0):
        f, axs  = plt.subplots(1)

        scattering_data = []
        for ray in self.hist["sca_angles"]:
            scattering_data.extend([a*RtoD for a in ray])

        w = np.ones_like(scattering_data) / len(scattering_data)
        hist, bins, _ = axs.hist(scattering_data, bins = 100, weights=w)

        axs.set_title("MS Scattering events angle histogram")


        #Make histogram of N random scattering angle given the relative cross section of Rayleigh and Mie scaterring.
        # N = 50000
        # f, axs  = plt.subplots(1)
        #
        # sca_angles_data = [self.GetScatteringAngle(cs_ray, cs_aer)*RtoD for i in range(N)]
        #
        # w = np.ones_like(sca_angles_data) / len(sca_angles_data)
        # hist, bins, _ = axs.hist(sca_angles_data, bins=self.atmosphere.profiles["sca_angle"]*RtoD, weights=w)
        #
        # phase_function = cs_ray * self.atmosphere.profiles["ray_Phase_Fct"] + cs_aer * self.atmosphere.profiles["aer_Phase_Fct"]
        # phase_function /= cs_ray + cs_aer
        # phase_function /= np.sum(phase_function)
        # # phase_function /= 2
        #
        # # print(np.sum(hist), np.sum(phase_function))
        #
        # # axs1 = axs.twinx()
        # axs.plot(self.atmosphere.profiles["sca_angle"]*RtoD, phase_function, "r")
        #
        # axs.set_title("MS scattering angle histogram")
