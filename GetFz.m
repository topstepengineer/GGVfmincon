function [Fz_fl, Fz_fr, Fz_rl, Fz_rr] = GetFz(setup, Fz0, FAero, accy, accx)

    dFz_lat = accy*(setup.mCar-setup.mUnsprung_f*2-setup.mUnsprung_r*2)*(setup.hCoG/setup.track);
    dFz_long = accx*(setup.mCar-setup.mUnsprung_f*2-setup.mUnsprung_r*2)*(setup.hCoG/setup.wheelbase);
    dFz_lat_f = setup.rMechbal*dFz_lat;
    dFz_lat_r = (1-setup.rMechbal)*dFz_lat;

    Fz_fl = Fz0(1) - dFz_lat_f - dFz_long + FAero(1)/2;
    Fz_fr = Fz0(2) + dFz_lat_f - dFz_long + FAero(1)/2;
    Fz_rl = Fz0(3) - dFz_lat_r + dFz_long + FAero(2)/2;
    Fz_rr = Fz0(4) + dFz_lat_r + dFz_long + FAero(2)/2;

end

