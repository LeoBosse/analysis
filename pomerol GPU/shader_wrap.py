
import os

import moderngl
import numpy as np
import imageio  # for output

import matplotlib.pyplot as plt

import moderngl as mgl
import struct

class ShaderWrap:

    def __init__(self, emission_origin, emission_map, uniforms, buffers, local_size=(1, 1, 1)):

        self.emission_origin = emission_origin
        self.emission_map = emission_map

        self.uniforms = uniforms
        self.buffer_data_list = buffers

        if self.emission_origin in ['grd', 'ground']:
            self.shader_file = 'grd_shader.glsl'
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



        self.LoadSource()

        self.CreateContext()
        self.CreateShader()

        self.WriteUniforms()
        self.InitBuffers()


    def LoadSource(self):
        ''' read and load the glsl code from file self.shader_file. '''
        with open(self.shader_file, 'r') as fp:
            self.src_code = fp.read()

        self.WriteConsts()

    def WriteConsts(self):
        for key, value in self.consts.items():
            self.src_code = self.src_code.replace(f"%%{key}%%", str(value))

        print("Writing uniforms lengths")
        for key, value in self.uniforms.items():
            if hasattr(value, '__iter__'): #if the value is a list of some sort, pass it as a texture
                self.src_code = self.src_code.replace(f"%%size_{key}%%", str(len(value)))
                # print(key, str(len(value)), f"%%size_{key}%%")

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
        # self.array_uniforms = []
        for key, value in self.uniforms.items():
            # if hasattr(value, '__iter__') and len(value) > 4: #if the value is a list of some sort, pass it as a texture
            #     value = np.array(value, dtype='float32')
            #     # self.array_uniforms.append(self.context.buffer(value))
            #     # print(key, value)
            #
            #     uniform = self.compute_shader.get(key, default='PROUT')
            #     self.compute_shader[key].write(data=value)
            #
            #     # print(key, value, self.compute_shader.get(key, default='PROUT').value)
            #
            #     # self.compute_shader.get(key, default=mgl.Uniform).value = value
            #
            # else:
            # print(key, value)
            uni = self.compute_shader.get(key, default=mgl.Uniform)
            uni.value = value
            # if uni != 'PROUT':
            #     # print('PROUT')
            # else:

            # print(key, value, self.compute_shader.get(key, default=mgl.Uniform).value)


    def InitBuffers(self):
        """Init the main input and output buffers to store the incoming and measured Stokes paramters. Iterable uniforms are created in WriteUniforms(). These buffers are 1D buffers the size of the input emission map to compute. the input V buffer is set with the radiance values of the emission map, while all others are set to zero."""
        # init buffers

        # test = np.zeros_like(self.emission_map.flatten(), dtype=np.float32)
        #
        # self.in_V_buffer_data    = np.array(self.emission_map.flatten(), dtype=np.float32)
        # self.in_Vcos_buffer_data = test #np.zeros_like(self.in_V_buffer_data, dtype=np.float32)
        # self.in_Vsin_buffer_data = test #np.zeros_like(self.in_V_buffer_data, dtype=np.float32)
        #
        # self.out_V_buffer_data    = test #np.zeros_like(self.in_V_buffer_data, dtype=np.float32)
        # self.out_Vcos_buffer_data = test #np.zeros_like(self.in_Vcos_buffer_data, dtype=np.float32)
        # self.out_Vsin_buffer_data = test #np.zeros_like(self.in_Vsin_buffer_data, dtype=np.float32)
        #
        #
        # self.in_V_buffer    = self.context.buffer(self.in_V_buffer_data)
        # self.in_Vcos_buffer = self.context.buffer(self.in_Vcos_buffer_data)
        # self.in_Vsin_buffer = self.context.buffer(self.in_Vsin_buffer_data)
        #
        # self.out_V_buffer    = self.context.buffer(self.out_V_buffer_data)
        # self.out_Vcos_buffer = self.context.buffer(self.out_Vcos_buffer_data)
        # self.out_Vsin_buffer = self.context.buffer(self.out_Vsin_buffer_data)

        self.buffer_list = []
        for b_data in self.buffer_data_list:
            # nb_lists = len(b_data)
            # data_len = len(b_data[0])

            # self.buffer_list.append(self.context.buffer(reserve= 4 * nb_lists * data_len, data=b_data))
            fmt_data = np.array(b_data, dtype=np.float32)
            # print(fmt_data, fmt_data.shape)
            self.buffer_list.append(self.context.buffer(fmt_data))

    def BindBuffers(self):
        """Bind the input and outpur stokes parameter buffers and the iterable uniforms."""
        #For Stokes paramters storage
        # self.in_V_buffer.bind_to_storage_buffer(0)
        # self.in_Vcos_buffer.bind_to_storage_buffer(1)
        # self.in_Vsin_buffer.bind_to_storage_buffer(2)
        # self.out_V_buffer.bind_to_storage_buffer(3)
        # self.out_Vcos_buffer.bind_to_storage_buffer(4)
        # self.out_Vsin_buffer.bind_to_storage_buffer(5)

        for i, buf in enumerate(self.buffer_list):
            buf.bind_to_storage_buffer(i)


    def Run(self):
        """Run the compute shader, Save the output buffers and release the context to allow the use of matplotlib."""

        self.BindBuffers()

        self.compute_shader.run(group_x=self.global_sizeX, group_y=self.global_sizeY, group_z=self.global_sizeZ)

        self.SaveBuffer()

        self.Release()


    def SaveBuffer(self):
        """Save the output buffers to useable numpy arrays."""

        print('RESULT IN')

        # test = np.frombuffer(self.in_V_buffer.read(), dtype=np.float32)
        # test = test.reshape(self.emission_map.shape)
        # print(test, test.shape)
        # test = np.frombuffer(self.in_Vcos_buffer.read(), dtype=np.float32)
        # test = test.reshape(self.emission_map.shape)
        # print(test, test.shape)
        # test = np.frombuffer(self.in_Vsin_buffer.read(), dtype=np.float32)
        # test = test.reshape(self.emission_map.shape)
        # print(test, test.shape)
        #
        # print('RESULT OUT')
        # test = np.frombuffer(self.out_V_buffer.read(), dtype=np.float32)
        # test = test.reshape(self.emission_map.shape)
        # print(test, test.shape)
        # test = np.frombuffer(self.out_Vcos_buffer.read(), dtype=np.float32)
        # test = test.reshape(self.emission_map.shape)
        # print(test, test.shape)
        # test = np.frombuffer(self.out_Vsin_buffer.read(), dtype=np.float32)
        # test = test.reshape(self.emission_map.shape)
        # print(test, test.shape)

        # for i in range(2):
        #     buffer = self.buffer_list[i]
        #     # print(buffer)
        #     data_shape = self.emission_map.shape
        #     test = np.frombuffer(buffer.read(), dtype=np.float32)
        #     # print(test, test.size)
        #     # print(test[::3].reshape(data_shape))
        #     # print(test[1::3].reshape(data_shape))
        #     # print(test[2::3].reshape(data_shape))
        #     test = test.reshape((data_shape[0], data_shape[1], 3))
        #     # print(test)


        self.output_shape = (self.global_sizeX, self.global_sizeY, self.local_sizeZ)

        out_V    = np.frombuffer(self.buffer_list[1].read(), dtype=np.float32)[0::3].reshape(self.output_shape)
        out_Vcos = np.frombuffer(self.buffer_list[1].read(), dtype=np.float32)[1::3].reshape(self.output_shape)
        out_Vsin = np.frombuffer(self.buffer_list[1].read(), dtype=np.float32)[2::3].reshape(self.output_shape)
        #
        # print(out_V)
        # print(out_Vcos)
        # print(out_Vsin)

        out_V    = np.sum(np.nan_to_num(out_V, nan=0), axis = 2)
        out_Vcos = np.sum(np.nan_to_num(out_Vcos, nan=0), axis = 2)
        out_Vsin = np.sum(np.nan_to_num(out_Vsin, nan=0), axis = 2)

        self.result = np.array([out_V, out_Vcos, out_Vsin])
        # print(self.result)
        # self.result = np.frombuffer(self.buffer_list[1].read(), dtype=np.float32).reshape((3, data_shape[0], data_shape[1]))
        # print("DEBUG results shape", self.result.shape)

    def Release(self):
        """Release the context. Mandatory to use matplotlib afterwards. Called automatically in the Run function."""
        self.context.release()
