import numpy
def windProfile(dx,dy,stations,inversion = 1000,topography = None):
    aa_c = 0.25
    ccc = 1.0
    eps = 1.e-10
    g = 9.81
    v_k = 0.4
    pi = numpy.pi
    tgr = 300.0
    aa2 = aa_c * (dx ** 2 + dy ** 2)
    stations["ust1"] = stations["velocity"](numpy.cos(stations["direction"]))
    stations["vst1"] = stations["velocity"](numpy.sin(stations["direction"]))
    stations["ust2"] = 1.5*stations["ust1"]
    stations["vst2"] = 1.5 * stations["vst1"]
    stations["temperature"] = tgr
    if topography is None:
        dhtop = 0
        dhmin = inversion
        sgmn = 0
        ixmin = 0
        iymin = 0

    beta = 0.2

    number_of_z_grid_points_below_0 = 3
    stab_ri_sign = 0

    Lc = h_canopy * (1 - l_p) / l_f
    ll = 2 * (beta ** 3) * Lc
    zz0 = ll / 0.4 * exp(-0.4 / beta)
    dd = ll / 0.4

    if (MetPar % is_const_wind) then
    U2h = MetPar % const_wind_sp
    theta = 270. - MetPar % const_wind_dir
    theta = theta * pi / 180.0
    else
    open(UNIT=29, file='h_stations.txt')
    REWIND(29)
    READ(29, *)
    h_stations
    close(29)

    sum_u = 0.0
    sum_v = 0.0
    total_sum = 0.0
    do
    i = 1, st % nst
    r2 = (x - st % xst(i)) ** 2 + (y - st % yst(i)) ** 2
    z_diff_sqr = (2 * h_canopy - h_stations(i)) ** 2
    weight_h = 1. / (1. + r2 / aa2)
    weight_v = 1. / (1. + ccc * z_diff_sqr / aa2)
    total_weight = weight_h * weight_v
    sum_u = sum_u + total_weight * st % ust(i, 1)
    sum_v = sum_v + total_weight * st % vst(i, 1)
    total_sum = total_sum + total_weight
    enddo
    uu = 1.0 * sum_u / (total_sum + norm_fac)
    vv = 1.0 * sum_v / (total_sum + norm_fac)
    U2h = sqrt((uu) ** 2 + (vv) ** 2)
    theta = atan2(vv, uu)
    endif
    Uh = (U2h * 0.4) / (beta * log((h_canopy + dd) / zz0))

    do
    iz = 1, wind % nwz
    z = max((iz - offset) * delta_z, 0.0)
    if (z.le.h_canopy) then
    Umag = Uh * exp(beta * (z - h_canopy) / ll)
    endif
    if (z.gt.h_canopy) then
    Umag = (1.0 * Uh * beta / 0.4) * log(1.0 * (z - h_canopy + dd) / zz0) ! in practice
    d = 0, so
    the
    altitude is above
    ground, not above
    zero
    plane !
    endif
    if (z.gt.2 * h_canopy) then
    Umag = U2h
    endif
    wind % u(ix, iy, iz) = cos(theta) * Umag
    wind % v(ix, iy, iz) = sin(theta) * Umag
    !if (h_canopy.gt.2) then
    !      print *, 'z = ', z, ' Umag = ', Umag, ' u = ', wind % u(ix, iy, iz), ' v = ', wind % v(ix, iy, iz)
    ! endif
    Umag = 0.0
    enddo
    enddo
    enddo
    InterpWind_2

    else
    !     ----interpolation
    on
    flow - surfaces - ----
    stab_ri_sign = 0
    InterpWind: do
    ix = 1, wind % nwx
    x = wind % dxru * (ix - 1)
    do
    iy = 1, wind % nwy
    y = wind % dyru * (iy - 1)
    h = InterpolateHight(x, y, topo)
    hx = DerivInterHight(1, 0, x, y, topo)
    hy = DerivInterHight(0, 1, x, y, topo)
    frc_htop = (h - topo % hmn) / dhtop
    delh = ht_max - h
    dis_lay = 0.5 * wind % dzru * delh
    do
    ist = 1, st % nst
    st % del_h(ist) = abs(h - st % hst(ist))
    enddo
    do
    iz = 1, wind % nwz
    pp = (iz - 1.5) * wind % dzru
    z = (ht_max - h) * pp
    z = max(z, 0.0)
    wei = 0.0
    weit = 0.
    D0 !!!!!.ef.added! ** *.EF.NISUYDIFF ** *.ef.diff
    wind % u(ix, iy, iz) = 0.0
    wind % v(ix, iy, iz) = 0.0
    wind % w(ix, iy, iz) = 0.0
    do
    i = 1, st % nst
    r2 = (x - st % xst(i)) ** 2 + (y - st % yst(i)) ** 2
    wei_h = 1. / (1. + r2 / aa2)  !horizontal
    weight
    wei_h = wei_h * st % west(i)
    zgal = max(z + h, z + st % hst(i))
    z_h = max(z, zgal - st % hst(i))
    !     determine
    the
    closest
    interpolation
    level
    if (i.le.st % nst_gr)then
    del_z = abs(zgal - st % hst(i) - 10.)
    ip_int = 1
    uu = st % ust(i, 1)
    vv = st % vst(i, 1)
    else
    if (iz.eq.1)then
    fct = log((zgal - h) / turb % z0(ix, iy)) / log(st % zprof(i, 2) / turb % z0(ix, iy))
    uu = st % ust(i, 2) * fct
    vv = st % vst(i, 2) * fct
    else
    do
    iprof = 1, st % nprof(i)
    zz = zgal - st % hst(i) - st % zprof(i, iprof)
    if (zz.lt.0.0)exit
    enddo
    if (zz.ge.0.0) then
    ip_int = st % nprof(i)
    uu = st % ust(i, ip_int)
    vv = st % vst(i, ip_int)
    else
    ip_int = iprof - 1
    del_z = zgal - st % hst(i) - st % zprof(i, ip_int)
    uu = st % ust(i, ip_int)
    vv = st % vst(i, ip_int)
    uu = st % ust(i, ip_int)
    vv = st % vst(i, ip_int)
    uu = uu + del_z / (st % zprof(i, ip_int + 1) - st % zprof(i, ip_int)) * &
    & (st % ust(i, ip_int + 1) - st % ust(i, ip_int))
    vv = vv + del_z / (st % zprof(i, ip_int + 1) - st % zprof(i, ip_int)) * &
    & (st % vst(i, ip_int + 1) - st % vst(i, ip_int))
    endif
    endif
    endif
    !phi = 1 / (sum: w_i)(sum: w_i * phi_i) =1 / a * b
    d_z = del_z - st % del_h(i) ** 3 / (st % del_h(i) ** 2 + z_h ** 2 + eps)
    wei_v = 1 / (1. + ccc * d_z ** 2 / aa2) !vertical
    weight
    wei = wei + wei_h * wei_v                                !! a = (sum: w_i)
    S2006, 10.1
    wind % u(ix, iy, iz) = wind % u(ix, iy, iz) + uu * wei_h * wei_v
    wind % v(ix, iy, iz) = wind % v(ix, iy, iz) + vv * wei_h * wei_v
    f_z = min((500.0 - z) / 300.0, 1.0)  !ls % factor
    for temperature
        weit = weit + wei_h * wei_v * f_z                          !! b = (sum: w_i * phi_i)
    MetPar % tmp(ix, iy, iz) = tt * wei_h * wei_v * f_z !!!!!!!!!  tt
    not defined !!???.ef.to
    check!!!!! to
    fix
    later - right
    now
    has
    no
    effect
    enddo
    wind % u(ix, iy, iz) = wind % u(ix, iy, iz) / wei
    wind % v(ix, iy, iz) = wind % v(ix, iy, iz) / wei
    MetPar % tmp(ix, iy, iz) = MetPar % tmp(ix, iy, iz) / weit
    !              print *, MetPar % tmp(ix, iy, iz), tt, wei_h, wei_v, f_z!.ef.debug
    wind % w(ix, iy, iz) = 0.0
    enddo
    enddo
    enddo
    InterpWind
    !  call
    wind_preplot(nisuy, topo, wind, t)

    # ifdef nodef
    !.ef.BUG !!!apperrs
    when
    running
    neutral!! CHECK
    this
    interpolation
    for STABLE stratification - introduces max height to climb according to brunt visala but has a bug !!!!
    !  if (MetPar % ind_stab.ne.1.and.stab_ri_sign.eq.1) then!.ef.diff was: if
        (MetPar % ind_stab.ne.1)
    then
    if (MetPar % ind_stab.ne.1) then
    StableLOOP: do
    ix = 1, wind % nwx
    do
    iy = 1, wind % nwy
    !     prepare
    parameters
    for stable layer
        !     calculates
        brunt - vassila
        frequency
        on
        grid
        points
    do
    iz = 1, wind % nwz
    MetPar % tmp(ix, iy, iz) = MetPar % t_gr(ix, iy) + MetPar % gam_gr * (iz - 1.5) * wind % dzru * delh
    MetPar % w_bv2(ix, iy, iz) = g * abs(MetPar % gam_gr) / MetPar % tmp(ix, iy, iz)
    enddo
    !     calculates
    possible
    vertical
    shift
    for each level
        do
        iz = 1, wind % nwz
    umag = sqrt(wind % u(ixmin, iymin, iz) ** 2 + wind % v(ixmin, iymin, iz) ** 2)
    tau_level = 1. / sqrt(MetPar % w_bv2(ix, iy, iz))
    fst(iz) = umag * tau_level
    devm = dhtop * (1 - (iz - 1.5) * wind % dzru)
    fst(iz) = min(fst(iz), devm)
    enddo
    do
    iz = 1, wind % nwz
    csi = (iz - 1.5) * wind % dzru
    siglev_h = ht_max * csi + (1 - csi) * h   !sigma
    level
    at
    point
    do
    izp = 1, wind % nwz
    sgs = sgmn(izp) + fst(izp) * frc_htop   !maximum
    height
    to
    climb
    if (sgs.gt.siglev_h) exit
    enddo
    if (sgs.le.siglev_h) izp=izp-1

    if (iz.gt.2)then
    !     interpolates
    fst
    sgsm1 = sgmn(izp - 1) + fst(izp - 1) * frc_htop
    rat = log(siglev_h / sgsm1) / log(sgs / sgsm1)
    fst_i = fst(izp - 1) + rat * (fst(izp) - fst(izp - 1))
    slp = fst_i / dhtop
    slp_rel = slp - (1 - (iz - 1.5) * wind % dzru)
    wind % w(ix, iy, iz) = (wind % u(ix, iy, iz) * hx + wind % v(ix, iy, iz) * hy) * slp_rel
    endif
    enddo
    wind % w(ix, iy, 2) = 0.0
    enddo
    enddo
    StableLOOP
    endif
    # endif

    !   Interpolation
    grid
    edges
    correction:
    for all stabilities
        EdgesZ: do
        ix = 1, wind % nwx
    do
    iy = 1, wind % nwy
    wind % u(ix, iy, 1) = 0.0 ! wind
    at
    z = 0 is zero
    wind % v(ix, iy, 1) = 0.0 ! wind
    at
    z = 0 is zero
    wind % w(ix, iy, 1) = 0.0 ! wind
    at
    z = 0 is zero
    wind % u(ix, iy, mz) = 2 * wind % u(ix, iy, mz - 1) - wind % u(ix, iy, mz - 2)
    wind % v(ix, iy, mz) = 2 * wind % v(ix, iy, mz - 1) - wind % v(ix, iy, mz - 2)
    wind % w(ix, iy, mz) = 2 * wind % w(ix, iy, mz - 1) - wind % w(ix, iy, mz - 2)
    enddo
    enddo
    EdgesZ
    EdgesY: do
    ix = 1, wind % nwx
    do
    iz = 1, wind % nwz
    wind % u(ix, 1, iz) = 2 * wind % u(ix, 2, iz) - wind % u(ix, 3, iz)
    wind % v(ix, 1, iz) = 2 * wind % v(ix, 2, iz) - wind % v(ix, 3, iz)
    wind % w(ix, 1, iz) = 2 * wind % w(ix, 2, iz) - wind % w(ix, 3, iz)
    wind % u(ix, mx, iz) = 2 * wind % u(ix, mx - 1, iz) - wind % u(ix, mx - 2, iz)
    wind % v(ix, mx, iz) = 2 * wind % v(ix, mx - 1, iz) - wind % v(ix, mx - 2, iz)
    wind % w(ix, mx, iz) = 2 * wind % w(ix, mx - 1, iz) - wind % w(ix, mx - 2, iz)
    enddo
    enddo
    EdgesY
    EdgesX: do
    iy = 1, wind % nwy
    do
    iz = 1, wind % nwz
    wind % u(1, iy, iz) = 2 * wind % u(2, iy, iz) - wind % u(3, iy, iz)
    wind % v(1, iy, iz) = 2 * wind % v(2, iy, iz) - wind % v(3, iy, iz)
    wind % w(1, iy, iz) = 2 * wind % w(2, iy, iz) - wind % w(3, iy, iz)
    wind % u(mx, iy, iz) = 2 * wind % u(mx - 1, iy, iz) - wind % u(mx - 2, iy, iz)
    wind % v(mx, iy, iz) = 2 * wind % v(mx - 1, iy, iz) - wind % v(mx - 2, iy, iz)
    wind % w(mx, iy, iz) = 2 * wind % w(mx - 1, iy, iz) - wind % w(mx - 2, iy, iz)
    enddo
    enddo
    EdgesX


    !     Calculate
    Ground
    FLUXES.ef. ! should
    be
    done
    on
    cartezian and not on
    masconst
    coords...
    UpdateFluxesGrnd: do
    ix = 1, wind % nwx
    x = wind % dxru * (ix - 1)
    do
    iy = 1, wind % nwy
    y = wind % dyru * (iy - 1)
    h = InterpolateHight(x, y, topo)
    dmh = ht_max - h
    dis_gr = 0.5 * wind % dzru * dmh
    cd = 0.4 / log(dis_gr / turb % z0(ix, iy)) ! sqrt
    of
    louis
    .13, Sqrt(Bulk
    Drag
    Coef)= v_star / u_mag = \kappa / ln[(z - z_d) / z0]

    !        print *, 'dzru=', wind % dzru, dmh, wind % nwz, dis_gr / exp(0.4 / cd)
    !        print *, 'cd=', cd, log(dis_gr / turb % z0(ix, iy)), dis_gr, dis_gr / turb % z0(ix, iy), ht_max, h
    u_mag = sqrt(wind % u(ix, iy, 2) ** 2 + wind % v(ix, iy, 2) ** 2)
    deltet = MetPar % gam_gr * dis_gr
    rib = 9.81 * dis_gr * deltet / (MetPar % t_gr(ix, iy) * u_mag ** 2) !  louis
    .11
    Ri_sign:
    if (rib.lt. - 0.005) then  !.ef.diff addition
    stab_ri_sign = 1
    else if (rib.gt.0.005) then!.ef.diff
    stab_ri_sign = -1
    else
    stab_ri_sign = 0
    endif
    Ri_sign
    !        if (rib.ge.0.0)then !.ef.diff
    if (stab_ri_sign.eq.1)then !.ef.diff
    fribm = 1. / (1. + 4.7 * rib) ** 2
    fribh = fribm
    else if (stab_ri_sign.eq. - 1)then
    cc = 9.4 * cd * sqrt(dis_gr / turb % z0(ix, iy))
    fribm = 1. - 9.4 * rib / (1 + cc * 7.4 * sqrt(abs(rib))) !louise
    .14
    fribh = 1. - 9.4 * rib / (1 + cc * 5.3 * sqrt(abs(rib))) !louise
    .15
    else !.ef.diff
    addition
    fribm = 1.
    fribh = 1.
    endif
    u_cof = max(u_mag, 1.5)  !!correction
    for low wind.ef.actually also for canopy, to check might be better with vstar threshhold !.ef..DEBUG.
    # ifdef min_turb_speed
    !! no
    practical
    effect, to
    check if can
    replace
    the
    correction
    for low wind
        turb % vstar(ix, iy) = max(u_cof * cd * sqrt(fribm), 0.2)  !.ef.minimal
        u_star, Hanna09:
        for urban(\sig_u_min=1 ) / 3=0.3
        # else
        turb % vstar(ix, iy) = u_cof * cd * sqrt(fribm)
        # endif
    if (MetPar % stability.eq."unstable") then
    turb % h_flux(ix, iy) = 0.1  ! heat
    flux, k
    m / s
    else
    turb % h_flux(ix, iy) = -cd ** 2 * u_mag * deltet * fribh / 0.74
    endif
    !print *, 'u_cof', u_cof, 'dis_gr', dis_gr, 'z0(ix,iy', turb % z0(ix, iy)
    !        if (stab_ri_sign.eq.0)turb % h_flux(ix, iy)=0.0 !neutral !.ef.diff
    enddo
    enddo
    UpdateFluxesGrnd
    print *, 'rib_sign=', stab_ri_sign
    return