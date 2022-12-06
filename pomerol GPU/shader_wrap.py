
import os

import moderngl
import numpy as np
import imageio  # for output

import matplotlib.pyplot as plt

import moderngl as mgl
import struct

class ShaderWrap:

    def __init__(self, emission_origin, emission_map, uniforms):

        self.emission_origin = emission_origin
        self.emission_map = emission_map

        self.uniforms = uniforms

        if self.emission_origin in ['grd', 'ground']:
            self.shader_file = 'grd_shader.glsl'
        else:
            raise('Shader type not recognize: use emission_origin=ground or sky.')

        self.global_sizeX = self.emission_map.shape[0] # Number of emission map pixels in distances
        self.global_sizeY = self.emission_map.shape[1] # Number of emission map pixels in azimuths
        self.global_sizeZ = 1  # Number of line of sight pixels in range
        # print(self.global_sizeX, self.global_sizeY, self.global_sizeZ)

        # Make sure these values are the same as in the shader .glsl code
        self.local_sizeX = 1
        self.local_sizeY = 1
        self.local_sizeZ = 1 #uniforms['los_nb_points']

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
        # print(self.src_code)

    def WriteUniforms(self):
        """Sets the uniform (or constant) parameters passed as a dictinionnary in <consts> """
        for key, value in self.uniforms.items():
            if hasattr(value, '__iter__'): #if the value is a list of some sort, pass it as a texture
                value = np.array([np.array(value).flatten()], dtype='float32')
                value = self.context.texture(size = value.shape, components = 1, data = value, dtype='f4')
                value.filter = (mgl.NEAREST, mgl.NEAREST)
                # value = text

            self.compute_shader.get(key, default=mgl.Uniform).value = value


    def CreateContext(self):
        """Create a openGL context to make the shader computations in. Don't forget to release it with context.release() once all computation are done and before using matplotlib."""
        self.context = mgl.create_standalone_context(require=430)

    def CreateShader(self):
        """Create a compute shader object from the context and the source code."""
        self.compute_shader = self.context.compute_shader(self.src_code)

    def InitBuffers(self):
        # init buffers
        test = np.zeros_like(self.emission_map.flatten(), dtype=np.float32)

        self.in_V_buffer_data    = np.array(self.emission_map.flatten(), dtype=np.float32)
        self.in_Vcos_buffer_data = test #np.zeros_like(self.in_V_buffer_data, dtype=np.float32)
        self.in_Vsin_buffer_data = test #np.zeros_like(self.in_V_buffer_data, dtype=np.float32)

        self.out_V_buffer_data    = test #np.zeros_like(self.in_V_buffer_data, dtype=np.float32)
        self.out_Vcos_buffer_data = test #np.zeros_like(self.in_Vcos_buffer_data, dtype=np.float32)
        self.out_Vsin_buffer_data = test #np.zeros_like(self.in_Vsin_buffer_data, dtype=np.float32)


        self.in_V_buffer    = self.context.buffer(self.in_V_buffer_data)
        self.in_Vcos_buffer = self.context.buffer(self.in_Vcos_buffer_data)
        self.in_Vsin_buffer = self.context.buffer(self.in_Vsin_buffer_data)

        self.out_V_buffer    = self.context.buffer(self.out_V_buffer_data)
        self.out_Vcos_buffer = self.context.buffer(self.out_Vcos_buffer_data)
        self.out_Vsin_buffer = self.context.buffer(self.out_Vsin_buffer_data)



    def Run(self):

        self.in_V_buffer.bind_to_storage_buffer(0)
        self.in_Vcos_buffer.bind_to_storage_buffer(1)
        self.in_Vsin_buffer.bind_to_storage_buffer(2)
        self.out_V_buffer.bind_to_storage_buffer(3)
        self.out_Vcos_buffer.bind_to_storage_buffer(4)
        self.out_Vsin_buffer.bind_to_storage_buffer(5)


        self.compute_shader.run(group_x=self.global_sizeX, group_y=self.global_sizeY, group_z=self.global_sizeZ)

        self.SaveBuffer()

        self.Release()


    def SaveBuffer(self):


        print('RESULT IN')

        test = np.frombuffer(self.in_V_buffer.read(), dtype=np.float32)
        test = test.reshape(self.emission_map.shape)
        print(test, test.shape)
        test = np.frombuffer(self.in_Vcos_buffer.read(), dtype=np.float32)
        test = test.reshape(self.emission_map.shape)
        print(test, test.shape)
        test = np.frombuffer(self.in_Vsin_buffer.read(), dtype=np.float32)
        test = test.reshape(self.emission_map.shape)
        print(test, test.shape)
        print('RESULT OUT')
        test = np.frombuffer(self.out_V_buffer.read(), dtype=np.float32)
        test = test.reshape(self.emission_map.shape)
        print(test, test.shape)
        test = np.frombuffer(self.out_Vcos_buffer.read(), dtype=np.float32)
        test = test.reshape(self.emission_map.shape)
        print(test, test.shape)
        test = np.frombuffer(self.out_Vsin_buffer.read(), dtype=np.float32)
        test = test.reshape(self.emission_map.shape)
        print(test, test.shape)


        out_V = np.frombuffer(self.out_V_buffer.read(), dtype=np.float32).reshape(self.emission_map.shape)
        out_Vcos = np.frombuffer(self.out_Vcos_buffer.read(), dtype=np.float32).reshape(self.emission_map.shape)
        out_Vsin = np.frombuffer(self.out_Vsin_buffer.read(), dtype=np.float32).reshape(self.emission_map.shape)

        self.result = np.array([out_V, out_Vcos, out_Vsin])
        # print(self.result, self.result.shape)

    def Release(self):
        self.context.release()


    # imgs = []
    # last_buffer = buffer_b
    # for i in range(FRAMES):
    #     toggle = True if i % 2 else False
    #     buffer_a.bind_to_storage_buffer(1 if toggle else 0)
    #     buffer_b.bind_to_storage_buffer(0 if toggle else 1)
    #
    #     # toggle 2 buffers as input and output
    #     last_buffer = buffer_a if toggle else buffer_b
    #
    #     # local invocation id x -> pixel x
    #     # work groupid x -> pixel y
    #     # eg) buffer[x, y] = gl_LocalInvocationID.x + gl_WorkGroupID.x * W
    #     compute_shader.run(group_x=H, group_y=1)
    #
    #     # print out
    #     output = np.frombuffer(last_buffer.read(), dtype=np.float32)
    #     output = output.reshape((H, W, 4))
    #     output = np.multiply(output, 255).astype(np.uint8)
    #     imgs.append(output)
    #
    # # if you don't want to use imageio, remove this section
    # out_path = f"{OUTPUT_DIRPATH}/debug.gif"
    # print("Writing GIF anim to", out_path)
    # imageio.mimwrite(out_path, imgs, "GIF", duration=0.15)
    #
