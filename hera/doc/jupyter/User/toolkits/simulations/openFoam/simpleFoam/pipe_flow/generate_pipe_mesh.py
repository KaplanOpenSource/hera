import os

from typing import Union
import classy_blocks as cb
import classy_blocks.construct.shapes.round as round_shapes

class pipe_mesher:
    @staticmethod
    def chop_round(shape : Union[round_shapes.RoundHollowShape, round_shapes.RoundSolidShape]):
        core_size = 0.5e-2
        bl_thickness =   0.7e-4
        shape.chop_axial(start_size=core_size)
        shape.chop_tangential(start_size=core_size)

        shape.chop_radial(end_size=bl_thickness, c2c_expansion=0.8)

    def create_pipe(self, target_dir):
        mesh = cb.Mesh() 

        pipe_shape  = cb.Cylinder(
            axis_point_1=[0, 0, 0], axis_point_2=[12, 0, 0], radius_point_1=[0, 4e-2, 0]
        )
        pipe_mesher.chop_round(pipe_shape)

        pipe_shape.set_start_patch("inlet")
        pipe_shape.set_end_patch("outlet")
        pipe_shape.set_outer_patch("walls")

        mesh.add(pipe_shape)

        mesh.modify_patch("walls", "wall")

        mesh.write(os.path.join(target_dir, "blockMeshDict"))