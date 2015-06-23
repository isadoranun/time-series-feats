def ReadLC_Catalina(lc):

    # Opening the lc
    data = pd.read_csv(lc,delim_whitespace=True,header=None)
    data.columns =['time','mag','error']

    data = data.sort(columns='time', ascending=True)

    mag = data.mag;
    time = data.time;
    error = data.error;


    return mag, time, error