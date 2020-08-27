interpolations
==============

Horizontal
**********
The interpolation is base on the IDW interpolation and use the elevation as additional weight function

.. code-block:: python

	interp(self, point, stations,topography=None, dx=20, dy=20, C=1000, D=5, Hsl=100, b=150)

Vertical
********

The vertical interpolations is based on the logaritmic and exponential profiles

.. code-block:: python

	windprofile(z, uref=3, href=24, he=24, lambdap=0.3, lambdaf=0.3,beta=0.3)
