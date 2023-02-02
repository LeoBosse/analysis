
import os

import numpy as np

import matplotlib.pyplot as plt

import moderngl as mgl

class ShaderWrap:

    def __init__(self, emission_origin, emission_map, local_size=(1, 1, 1)):

        self.emission_origin = emission_origin
        self.emission_map = emission_map
        #
        # self.uniforms = uniforms
        # self.buffer_data_list = buffers

        if self.emission_origin in ['grd', 'ground']:
            self.shader_file = 'grd_shader.glsl'
        elif self.emission_origin in ['sky']:
            self.shader_file = 'sky_shader.glsl'
        else:
            raise('Shader type not recognize: use emission_origin=ground or sky.')

        self.global_sizeX = self.emission_map.shape[0] # Number of emission map pixels in distances
        self.global_sizeY = self.emission_map.shape[1] # Number of emission map pixels in azimuths
        self.global_sizeZ = 1  # Number of line of sight pixels in range
        # print(self.global_sizeX, self.global_sizeY, self.global_sizeZ)

        # Make sure these values are the same as in the shader .glsl code
        self.local_sizeX = local_size[0]
        self.local_sizeY = local_size[1]
        self.local_sizeZ = local_size[2]

        self.consts = {'local_sizeX': self.local_sizeX,
                       'local_sizeY': self.local_sizeY,
                       'local_sizeZ': self.local_sizeZ
        }

        self.buffer_list = []

        self.CreateContext()
        self.LoadSource()
        self.CreateShader()


    def Prepare(self, uniforms, in_buffers, out_buffers):

        self.uniforms = uniforms
        
        self.nb_in_buffers = len(in_buffers)
        self.buffer_data_list = in_buffers
        self.buffer_data_list.extend(out_buffers)

        self.WriteUniforms()

        if self.buffer_list == []:
            self.InitBuffers()
            self.BindBuffers()
        else:
            self.WriteBuffers()
            # self.BindBuffers()


    def LoadSource(self):
        ''' read and load the glsl code from file self.shader_file.'''
        with open(self.shader_file, 'r') as fp:
            self.src_code = fp.read()

        self.WriteConsts()

    def WriteConsts(self):
        for key, value in self.consts.items():
            self.src_code = self.src_code.replace(f"%%{key}%%", str(value))

        # print(self.src_code)

    def CreateContext(self):
        """Create a openGL context to make the shader computations in. Don't forget to release it with context.release() once all computation are done and before using matplotlib."""
        self.context = mgl.create_standalone_context(require=440)

        # print("DEBUG self.context.info", self.context.info)
        # print("DEBUG self.context.info", self.context.info['GL_MAX_COMPUTE_UNIFORM_COMPONENTS'])
        # print("DEBUG self.context.info", self.context.info['GL_MAX_UNIFORM_LOCATIONS'])

    def CreateShader(self):
        """Create a compute shader object from the context and the source code."""
        self.compute_shader = self.context.compute_shader(self.src_code)

        # print("All uniforms")
        # for u in self.compute_shader:
        #     print(u)

    def WriteUniforms(self):
        """Sets the uniform parameters as in self.uniforms."""
        for key, value in self.uniforms.items():
            self.compute_shader[key].value = value

        # print('-------------shader uniforms-------------')
        # for key in self.compute_shader:
        #     print(key, self.compute_shader[key].value)

    def InitBuffers(self):
        """Init the main input and output buffers to store the incoming and measured parameters."""
        self.buffer_list = []
        for b_data in self.buffer_data_list:
            fmt_data = np.array(b_data, dtype=np.float32)
            self.buffer_list.append(self.context.buffer(fmt_data))


    def WriteBuffers(self):
        for i, b_data in enumerate(self.buffer_data_list):
            fmt_data = np.array(b_data, dtype=np.float32)
            # self.buffer_list[i].clear()
            self.buffer_list[i].write(fmt_data)


    def BindBuffers(self):
        """Bind the input and output buffers to a single int identifier. The buffer_list MUST be the same order than their declaration in the shader source code."""

        for i, buf in enumerate(self.buffer_list):
            buf.bind_to_storage_buffer(i)


    def Run(self):
        """Run the compute shader, Save the output buffers."""

        self.compute_shader.run(group_x=self.global_sizeX, group_y=self.global_sizeY, group_z=self.global_sizeZ)

        self.SaveBuffer()
        # self.OrphanBuffers()

        # self.Release()


    def SaveBuffer(self, buffer_ID = None):
        """Save the output buffers to useable numpy arrays."""


        if buffer_ID is None:
            buffer_ID = self.nb_in_buffers

        self.output_shape = (self.global_sizeX, self.global_sizeY, self.local_sizeZ)

        out_V    = np.frombuffer(self.buffer_list[buffer_ID].read(), dtype=np.float32)[0::3].reshape(self.output_shape)
        out_Vcos = np.frombuffer(self.buffer_list[buffer_ID].read(), dtype=np.float32)[1::3].reshape(self.output_shape)
        out_Vsin = np.frombuffer(self.buffer_list[buffer_ID].read(), dtype=np.float32)[2::3].reshape(self.output_shape)


        # print('-------------RESULT 000-------------')
        # print(out_V)
        # print(out_Vcos)
        # print(out_Vsin)

        out_V    = np.sum(np.nan_to_num(out_V,    nan=0), axis = 2)
        out_Vcos = np.sum(np.nan_to_num(out_Vcos, nan=0), axis = 2)
        out_Vsin = np.sum(np.nan_to_num(out_Vsin, nan=0), axis = 2)

        self.result = np.array([out_V, out_Vcos, out_Vsin])
        # print(self.result)
        # print(self.result)
        # self.result = np.frombuffer(self.buffer_list[1].read(), dtype=np.float32).reshape((3, data_shape[0], data_shape[1]))
        # print("DEBUG results shape", self.result.shape)

    def Release(self):
        """Release the context. Mandatory to use matplotlib afterwards."""
        self.context.release()
    # def OrphanBuffers(self):
    #     for i in range(len(self.buffer_list)):
    #         self.buffer_list[i].orphan()

        del self.context



class ShaderWrapMS(ShaderWrap):

    def __init__(self, N_rays, sca_limit):
        self.shader_file = 'MS_shader.glsl'

        self.size =  int(N_rays ** (1/6.)) + 1

        self.N_rays = self.size ** 6

        self.global_sizeX = self.size # Number of emission map pixels in distances
        self.global_sizeY = self.size # Number of emission map pixels in azimuths
        self.global_sizeZ = self.size  # Number of line of sight pixels in range
        # print(self.global_sizeX, self.global_sizeY, self.global_sizeZ)

        # Make sure these values are the same as in the shader .glsl code
        self.local_sizeX = self.size
        self.local_sizeY = self.size
        self.local_sizeZ = self.size

        self.consts = {'local_sizeX': self.local_sizeX,
                       'local_sizeY': self.local_sizeY,
                       'local_sizeZ': self.local_sizeZ,
                       'scattering_limit': sca_limit
        }

        self.buffer_list = []

        self.CreateContext()
        self.LoadSource()
        self.CreateShader()


    def SaveBuffer(self, buffer_ID = None):
        """Save the output buffers to useable numpy arrays."""
        if buffer_ID is None:
            buffer_ID = self.nb_in_buffers

        out_V    = np.frombuffer(self.buffer_list[buffer_ID].read(), dtype=np.float32)[0::3]
        out_Vcos = np.frombuffer(self.buffer_list[buffer_ID].read(), dtype=np.float32)[1::3]
        out_Vsin = np.frombuffer(self.buffer_list[buffer_ID].read(), dtype=np.float32)[2::3]

        # print('-------------RESULT 000-------------')
        # print(out_V)
        # print(out_Vcos)
        # print(out_Vsin)

        sum_out_V    = np.sum(np.nan_to_num(out_V,    nan=0))
        sum_out_Vcos = np.sum(np.nan_to_num(out_Vcos, nan=0))
        sum_out_Vsin = np.sum(np.nan_to_num(out_Vsin, nan=0))

        self.result = np.array([out_V, out_Vcos, out_Vsin])
        # print(self.result)
