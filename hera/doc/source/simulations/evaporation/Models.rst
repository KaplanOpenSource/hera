Evaporation Models
==================

Evaporation Rate
----------------

US
..

This model uses the next equations for the exaporation rate:

.. math::
    J=\frac{K_m M_w P_s (T)}{RT}=flux\;rate\;\left(\frac{g}{m^2s}\right) \\
    K_m=
    \begin{cases}
    0.664Re^{-0.5}Sc^{-2/3}u\;\;\;Re\leq20,000\\
    0.0336Re^{-0.2}Sc^{-2/3}u\;\;\;Re>20,000
    \end{cases}\\
    M_w=molecular\;weight\;\left(\frac{g}{mol}\right)\\
    P_s=vapor\;pressure\left\;(atm\right)\\
    R=gas\;constant\\
    T=temperature\;(K)\\
    Re=\frac{lu}{\nu_air}=Reynolds\;number\\
    Sc =\frac{\nu_air}{D}=Schmidt\;number\\
    u = wind\;velocity\;\left(\frac{m}{s}\right)\\
    l=pool\;length\;downwind\;(m)\\
    \nu_air = air\;kinematic\;viscocity\;\left(\frac{m^2}{s}\right)\\
    D = diffusivity\;\left(\frac{m^2}{s}\right)\\

Diffusion
---------

FSG
...

This model uses the next equations for the molecular diffusion:

.. math::
    D =10^{-3}\frac{T^{1.75}\sqrt{1/M_air+1/M_agent}}{\left(V_air^{1/3}+V_agent^{1/3}\right)^2}
    =diffusivity\;\left(\frac{m^2}{s}\right)\\
    T=temperature\;(K)\\
    M=molecular\;weight\;\left(\frac{g}{mol}\right)\\
    V=liquid\;molar\;volume\;\left(\frac{ml}{mol}\right)\\

EPA
...

This model uses the next equations for the molecular diffusion:

.. math::
    D =0.0000409\frac{T^{1.9}\sqrt{1/M_air+1/M_agent}}{M_agent^{1/3}}
    =diffusivity\;\left(\frac{m^2}{s}\right)\\
    T=temperature\;(K)\\
    M=molecular\;weight\;\left(\frac{g}{mol}\right)\\

Dynamic viscocity
-----------------

The kinematic viscocity of the air is calculated using its density and
dynamic viscocity.
The density is calculated using the perfect gas law.
The dynamic viscocity is calculated using one of the next models:

Power law
.........

.. math::
    \mu =1.8205\cdot10^{-5} \sqrt{\frac{T}{293}}
    =dynamic\;viscocity\;\left(\frac{kg}{ms}\right)\\
    T=temperature\;(K)\\

Vapor Pressure
--------------
Antoine equation
................
The vapor pressure may be calculated using Antoine equation.

.. math::
    log(P) =A-\frac{B}{T-C}\\
    T=temperature\;(K)\\

A, B, C are constants that should be specified in the agent's physical properties dictionary.
The units of the pressure P should be specified there too.
