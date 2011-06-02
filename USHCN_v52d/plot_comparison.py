
import numpy as np
import pylab

def read_comp(station_id, src):
    '''
    
    '''
    
    if src == "hom":
        fn = "src_codes/benchmark/monthly/WMs.52d/%s_avg.WMs.52d" % station_id
        year_bit = 0
        month_indices = range(1, 13)
    elif src == "raw":
        fn = "src_codes/benchmark/monthly/raw/%s_avg.raw" % station_id
        year_bit = 1
        month_indices = range(2,14)

    print "...Reading %s" % fn
    
    datafile = open(fn, "r")    
    station_data = dict()

    years = []
    datas = []
    for l in datafile:
        print l.strip()
        bits = l.strip().split()

        year = int(bits[year_bit][-4:])
        years.append(year)
        print "...", year

        monthly_bits = np.array(bits)[month_indices]
        try:
            data_from_file = 0.1*np.array(map(int, monthly_bits))
        except (ValueError, ):
            print "......Caught missing data - replacing with null values"
            monthly_data = []
            
            for value in monthly_bits:
                if value.isdigit():
                    monthly_data.append(0.1*int(value))
                else:
                    monthly_data.append(None)
            data_from_file = np.ma.array(monthly_data)
            
        datas.append(data_from_file)
        print "...", data_from_file, "\n"
        
    years = np.ma.array(years)
    datas = np.ma.array(datas, dtype='float')

    station_data['year'] = years
    station_data['data'] = datas

    all_data_monthly = np.concatenate([datas[i,:] for i in range(len(years))])
    all_data_monthly = np.ma.array(all_data_monthly, dtype='float')
    
    station_data['monthly'] = all_data_monthly

    return station_data             
             
if __name__ == "__main__":

    station_ids = range(830000, 830008)
    all_station_data = dict()
    for station_id in station_ids:
        raw_station_data = read_comp(station_id, "raw")
        hom_station_data = read_comp(station_id, "hom")
 
        all_station_data[station_id] = { 'raw' : raw_station_data['monthly'],
                                         'hom' : hom_station_data['monthly'], }
        
    numplots = len(station_ids)
    pylab.figure(figsize=[12*2, numplots*2])
    pylab.subplots_adjust(hspace=.6)
    params = { 'legend.fontsize': 8,
               'xlabel.fontsize': 8,
               'ylabel.fontsize': 8, }
    pylab.rcParams.update(params)              
    
    for p in range(numplots):
        
        pylab.subplot(numplots/2, 2, p+1)
        station_id = station_ids[p]
        print p, station_id
        station_data = all_station_data[station_id]
        raw = station_data['raw']
        hom = station_data['hom']
        
        pylab.plot(raw)
        pylab.plot(hom)

        pylab.title("Station - %6d" % station_id)
        pylab.xticks( ([0, 240, 480, 720, 960, 1199]), map(str, [1900, 1920, 1940, 1960, 1980, 2000]),
                      rotation=17 )        
        
        adjustments = []
        for (hom_datum, raw_datum) in zip(hom, raw):
            if hom_datum and raw_datum:
                adjustment = hom_datum-raw_datum
            else:
                adjustment = None
            adjustments.append(adjustment)
        adjustments = np.ma.array(adjustments)
        pylab.plot(adjustments)
               
        if p == 0:
            pylab.legend( ('raw', 'homogenized', 'adjustments') )        
    pylab.show()

