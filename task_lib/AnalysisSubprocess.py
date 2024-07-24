from XATOMS.postprocessors.ovito_calculators import OvitoCalculators
from XATOMS.postprocessors import compute_twist, compute_diffraction
import math

def AnalysisSubprocess(comm, input_params):

        nprocs = comm.Get_size()
        rank = comm.Get_rank()
        MAX_STRIDES = math.ceil(input_params["total_number_of_timesteps"]/input_params["i_o_freq"]/(nprocs - input_params['md_procs']))
        stride = 0 #

        while stride<MAX_STRIDES:
                
       # receive LAMMPS Snapshot


                timestep = (stride*(nprocs-input_params["md_procs"]) + (rank-input_params["md_procs"]))*input_params["i_o_freq"]

                # print (f'currently at rank {rank} and stride {stride} at time step: {timestep} BEFORE receving')
                
                if timestep<=input_params["total_number_of_timesteps"]:
                        lmp_snapshot = comm.recv(source=0)
                        # print (f'currently at rank {rank} and stride {stride} AFTER receving')
                else:
                        return
                
                # Run analysis for each snapshot
                
                data = OvitoCalculators(lmp_snapshot=lmp_snapshot)


                if 'dxa_analysis' in input_params.keys():

                        if input_params["dxa_analysis"]:
                                data.calculate_dxa(export_formats=input_params["dxa_output_formats"], dxa_line_sep=input_params["dxa_line_sep"], lattice="lattice")
                
                if input_params["full_trajectory_dump"]:
                        data.dump_trajectories(export_formats=input_params["trajectory_output_format"])
                
                if 'compute_twist' in input_params.keys():
                         
                        twist_angle, err1, err2 = compute_twist.get_interlayer_twist(input_params['compute_twist']['cutoff'], \
                                                                                     input_params['compute_twist']['id_1'], input_params['compute_twist']['id_2'], \
                                                                                     input_params['compute_twist']['reference_particle_type'], input_params['compute_twist']['num_iter'])
                        with open(f'twist_{lmp_snapshot.timestep}', 'w') as file:
                                file.write(f'time-step twist-angle fit_err_layer_1 fit_err_layer_2\n')
                                file.write(f'{lmp_snapshot.timestep} {twist_angle} {err1} {err2}')
        
                
                if 'compute_Laue_Diffration' in input_params.keys():
                        
                        if input_params['compute_Laue_Diffraction']:
                                filetag = str(lmp_snapshot.timestep)
                                compute_diffraction.get_laue_pattern(data, filetag)

                if 'compute_xrd' in input_params.keys():

                        if input_params['compute_xrd']:
                                filetag = str(lmp_snapshot.timestep)
                                compute_diffraction.get_xrd_pattern(data, filetag)

                        
                stride+=1


        return