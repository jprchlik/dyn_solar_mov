import matplotlib
#fixes multiprocess issue
matplotlib.use('agg')
import sunpy.map

from sunpy.cm import cm
from sunpy.net.helioviewer import HelioviewerClient
from sunpy.physics import solar_rotation
from sunpy.net import hek

import subprocess
import glob
import os
import sys, getopt
import numpy as np
from multiprocessing import Pool

from datetime import date,datetime
from datetime import timedelta as dt

import matplotlib.pyplot as plt
import astropy.units as u

hv = HelioviewerClient()

#create variables with defaults for my machine
#starting directory
stard = '/Volumes/Pegasus/jprchlik/video_wall_testing/'
sdir = stard
#add information about named AR
addinfo = True
#dimensions of the image
#Initially for videowall in Solar B
w0 = 1900
h0 = 1144
#movie frame rate
frate = 25
#output movie file
outmov = 'dailysun'
#hour between aia update
dh = 0.
#hour between aia update
dm = 40.
#set dpi
dpi = 300
#number of processors for downloading
nproc = 10
#information on how long to run movie
span = 2. #days to run the movie over
cadence = 6. #in minutes
minweek = 10080. #minutes in a week
#day array
dayarraya = []
#list of objects
objelist = []
#dict of flares
eventobj = {}


def get_file(ind):
    global hv,dayarraya
    global stard,addinfo,w0,h0,frate,outmov,dh,dm,dpi,nproc,span,cadence,minweek,sdir
#Source ID 11 is AIA 193
#    fils = ind.replace(':','_').replace('/','_').replace(' ','_')+'_AIA_193.jp2'
    try:
        meta = hv.get_closest_image(ind,sourceId=11)
    except:
        print 'Failed to grab '+ind
        return
    date = meta['date'].strftime('%Y_%m_%d__%H_%M_%S')
    mils = float(meta['date'].strftime('%f'))/1000.
    fstr = date+'_*__SDO_AIA_AIA_193.jp2'.format(mils).replace(' ','0')
#    test = os.path.isfile(sdir+'/raw/'+fstr)
# check to see if the working directory file already exists because I will remove the raw file
    test = glob.glob(sdir+'/working/'+fstr.replace('jp2','png'))
    traw = glob.glob(sdir+'/raw/'+fstr)
#test to make sure cadence matches 
#Does not work as of (2016/09/26 J. Prchlik)
#set up to subtract from the larger because datetime is funny if you subtract the past from the future
#date time object of asking date
    dtind = datetime.strptime(ind,'%Y/%m/%d %H:%M:%S')
    difft = abs(dtind-meta['date']).seconds #difference in time in seconds
# set the tolerance for the cadence to be a valid image
#currently setting to be a static 30s
    testcad = difft < cadence*60./4. 
#    testcad = True
        
    if ((len(test) == 0) & (len(traw) == 0) & (testcad)):
#    if len(test) == 0:
#        filep = hv.download_png(ind,0.3,"[11,1,100]",directory=sdir+'/raw',x0=0,y0=0,width=int(rv*8192),height=8192,watermark=False)#,y1=-1200,y2=1200,x1=-1900,x2=1900)
        try:
            filep = hv.download_jp2(ind,sourceId="11",directory=sdir+'/raw',clobber=True)#,y1=-1200,y2=1200,x1=-1900,x2=1900)
            check, img = qual_check(filep)
            if check != True:
    #do not create png file if the image does not pass quality checks
                print filep+' did not pass quality checks'
                os.remove(filep)
        except:
            print 'Failed to grab '+ind
            return
#for j,i in enumerate(dayarray):
def qual_check(filep):
#read JPEG2000 file into sunpymap
    img = sunpy.map.Map(filep)
#Level0 quality flag equals 0 (0 means no issues)
    lev0 = img.meta['quallev0'] == 0
#check level1 bitwise keywords (http://jsoc.stanford.edu/doc/keywords/AIA/AIA02840_K_AIA-SDO_FITS_Keyword_Document.pdf)
    lev1 = np.binary_repr(img.meta['quality']) == '1000000000000000000000000000000'
#check to see if it is a calibration image
#This keyword changed after AIA failure
#    calb = np.binary_repr(img.meta['aiftsid']) == '1010000000000000'
#check that both levels pass and it is not a calibration file
    check = ((lev0) & (lev1))# & (calb)) 

    return check,img

#return image extent for giving coordinates to image
def img_extent(img):
# get the image coordinates in pixels
    px0 = img.meta['crpix1']
    py0 = img.meta['crpix2']
# get the image coordinates in arcsec 
    ax0 = img.meta['crval1']
    ay0 = img.meta['crval2']
# get the image scale in arcsec 
    axd = img.meta['cdelt1']
    ayd = img.meta['cdelt2']
#get the number of pixels
    tx,ty = img.data.shape
#get the max and min x and y values
    minx,maxx = px0-tx,tx-px0
    miny,maxy = py0-ty,ty-py0
#convert to arcsec
    maxx,minx = maxx*axd,minx*axd
    maxy,miny = maxy*ayd,miny*ayd

    return maxx,minx,maxy,miny



#reformat file to be in 1900x1200 array and contain timetext
def format_img(i):
    global stard,addinfo,w0,h0,frate,outmov,dh,dm,dpi,nproc,span,cadence,minweek,sdir
    global eventobj,objelist,addinfo,dayarraya
    global cadence
    filep = dayarraya[i]
# test to see of symbolic link exits for given number and replace it if it does
    symli = 'working/symlinks/seq{0:4d}.png'.format(i).replace(' ','0')
#    testsym = os.path.isfile(symli)
#    if testsym:
#        os.remove(testsym)	

    #scale up the images with increasing dpi
    sc = dpi/100
#output file same name as raw file except png
    outfi = filep.replace('raw','working').replace('jp2','png')
    test = os.path.isfile(outfi)

#test to see if working file exists
    if test == False:
#  run quality checks on the data
        check,img = qual_check(filep)
#        if check:
        print 'Modifying file '+filep
#        img = sunpy.map.Map(filep)
        fig,ax = plt.subplots(figsize=(sc*float(w0)/float(dpi),sc*float(h0)/float(dpi)))
        fig.set_dpi(dpi)
        fig.subplots_adjust(left=0,bottom=0,right=1,top=1)
        ax.set_axis_off()
#return extent of image
        maxx,minx,maxy,miny = img_extent(img)

#plot the image in matplotlib
        ax.imshow(img.data,interpolation='none',cmap=cm.sdoaia193,vmin=0,vmax=255,origin='lower',extent=[minx,maxx,miny,maxy])
        ax.set_axis_bgcolor('black')
        ax.text(-2000,-1100,'AIA 193 - '+img.date.strftime('%Y/%m/%d - %H:%M:%S')+'Z',color='white',fontsize=36,zorder=50,fontweight='bold')
        tstart = img.date.strftime('%Y/%m/%d %H:%M:%S')
#list of indices to remove
#        remove_list = []
#loop over event dictionary at start of process
        for j,i in enumerate(objelist):
            arname = i['ar_noaanum']
#start time of events
            sevent = datetime.strptime(i['event_starttime'],'%Y-%m-%dT%H:%M:%S')
#get ar coordiantes
            x,y = i['hpc_coord'].replace('POINT(','').replace(')','').split(' ')
            x,y = float(x),float(y)
            x,y = solar_rotation.rot_hpc(x*u.arcsec,y*u.arcsec,sevent,img.date)
            x,y = x.value,y.value
#previous image x,y value
#            xp,yp = solar_rotation.rot_hpc(x*u.arcsec,y*u.arcsec,sevent,img.date+dt(minutes=-cadence*2))
#            xp,yp = xp.value,yp.value
#Only plot AR if is on the front side of the sun keep 1 cadence past end time
#            if img.date <= eventobj[str(arname)+'_endtime']+dt(minutes=30.*cadence):
#Three hours is a buffer for HEK to update with the next day's observations
            if img.date <= eventobj[str(arname)+'_endtime']+dt(hours=3):
                ax.text(x,y,str(j+1),alpha=0.45,zorder=150,color='black',fontsize=18)
                ax.text(-2000,1200-(j+1)*200.,'{0:1d}. AR {1:5d}'.format(j+1,arname),color='white',fontsize=26)
                ax.text(-2000,1200-(j+1)*200.-50.,'{1}'.format(i['ar_numspots'],eventobj[str(arname)+'_sflare']),color='white',fontsize=16)
                ax.text(-2000,1200-(j+1)*200.-100.,'Sunspots = {0:3d}'.format(i['ar_numspots'],eventobj[str(arname)+'_sflare']),color='white',fontsize=16)
#            else:
    #add index to removal list
   #                 remove_list.append(j)
    #causes flashing do no do
    #remove object when it goes over the limb
    #        for j in remove_list:
    #            try:
    #                objelist = np.delete(objelist,j)
    #            except:
    #                print 'Event Already removed'
    
    
            fig.savefig(outfi,edgecolor='black',facecolor='black',dpi=dpi)
            plt.clf()
            plt.close()
#remove previous raw dir file 
        os.remove(filep)
# no matter what change the symbolic link number
    os.symlink(outfi,symli)
    return

	
#remove symbolic links in parallel
def remove_sym(files):
    global stard,addinfo,w0,h0,frate,outmov,dh,dm,dpi,nproc,span,cadence,minweek,sdir
    global dayarraya
    os.remove(files)
    return
#remove all previous symbolic links

def main(argv):
    global stard,addinfo,w0,h0,frate,outmov,dh,dm,dpi,nproc,span,cadence,minweek,sdir
    global dayarraya,objelist,eventobj
    #arugments the prgram to read
    inargs1 = 'hf:i:w:h:r:m:lh:lm:d:p:s:c:'
    snargs1 = inargs1[1:].split(':')
    inargs2 = ["sdir","addinfo","width","height","frate","outmov","laghour","lagmin","dpi","nproc","span","caden"]
    helpinfo = "run_video_wall.py is a program that creates a movie using solar AIA 193\n"
    helpinfo = helpinfo+"The command takes many arguments but many values are filled in by default if you are running for the video wall at SAO\n"
#    helpinfo = helpinfo+"python run_video_wall.py -d <sdir> -i <addinfo> -w <width> -h <height> -r <frate> -m <outmov> -lm <laghout>"
    helpinfo = helpinfo+"python run_video_wall.py"
#basic help creation
    for i in range(len(inargs2)):
        helpinfo = helpinfo+' -'+snargs1[i]+' <--'+inargs2[i]+'>'
    helpinfo = helpinfo+'\n'

    

#Descripted information about each keyword works
    argsdes = ["The directory where all of the file creation occurs. This directory must already exist, but the program will create all subdirectories",
	           "Adds information about AR regions including position, number of flares, and number of sunspots (Default = True)",
               "Width of the image and created movie in pixels (Default ={0:5d})".format(w0),
               "Height of the image and created movie in pixels (Default ={0:5d})".format(h0),
               "Frame rate of the output movie (Default ={0:2d})".format(frate),
               "Output movie filename (do not add extension) (Default = {0})".format(outmov),
               "Lag time between the current UTC time and the output movie time in UTC in hours (Default ={0:5.2f})".format(dh),
               "Lag time between the current UTC time and the output movie time in UTC in minutes (Default ={0:5.2f})".format(dm),
               "Dots per inch in the created images (Default ={0:3d})".format(dpi),
               "Number of processors to used (Default ={0:3d})".format(nproc),
               "The lengths of the movie in days (Default ={0:2.0f})".format(span),
               "Cadence of the movie in minutes (Default ={0:5.2f})".format(cadence)]

    for i in range(len(inargs2)):
        helpinfo = helpinfo+' -'+snargs1[i]+' <--'+inargs2[i]+'> : '+argsdes[i]+'\n'

    try:
        opts, args = getopt.getopt(argv,inargs1,inargs2)
    except getop.GetoptError:
        print helpinfo
        sys.exit(2)
    inargs1 = 'hf:i:w:h:r:m:lh:lm:d:p:s:c:'
    snargs1 = inargs1[1:].split(':')
    inargs2 = ["sdir","addinfo","width","height","frate","outmov","laghour","lagmin","dpi","nproc","span","caden"]

    #load in user values
    for opt, arg in opts:
        if opt == '-h':
            print helpinfo
            sys.exit()
        elif opt in ("-f","--sdir"):
            sdir  = str(arg)
            stard = str(arg)
        elif opt in ("-i","--addinfo"):
            addinfo = arg
            if addinfo.upper() == 'FALSE':
                addinfo = False
            else:
                addinfo = True
        elif opt in ("-w","--width"):
            w0 = int(arg)
        elif opt in ("-h","--height"):
            h0 = int(arg)
        elif opt in ("-r","--frate"):
            frate = int(arg)
        elif opt in ("-m","--outmov"):
            outmov = arg
        elif opt in ("-lh","--laghour"):
            dh = float(arg)
        elif opt in ("-lm","--lagmin"):
            dm = float(arg)
        elif opt in ("-d","--dpi"):
            dpi = int(arg)
        elif opt in ("-p","--nproc"):
            nproc = int(arg)
        elif opt in ("-s","--span"):
            span = float(arg)
        elif opt in ("-c","--caden"):
            cadence = float(arg)
            
    
    
    #video wall ratio
    rv = float(w0)/float(h0)
    
    #helioviewer client
    
    
    
    
    
    #output image dimensions 
    x1 = -w0
    x2 =  w0
    y1 = -h0
    y2 =  h0
    
    
    #have last time be now it UTC
    now = datetime.utcnow()
    eday = now
    
    #change the directory to the main directory
    os.chdir(stard)
    #The directory which will contain the raw png files
    sdir = stard
    try:
    #    os.mkdir(sdir)
        os.mkdir(sdir+'/raw')
        os.mkdir(sdir+'/working')
        os.mkdir(sdir+'/working/symlinks')	
        os.mkdir(sdir+'/final')
    except OSError:
    #remove symbolic links from before
        symfiles = glob.glob(sdir+'/working/symlinks/*png')
        poolsym = Pool(processes=nproc)
        poolsym.map(remove_sym,symfiles)
        poolsym.close()
    #remove leftover raw files in case of error in previous run
        rawfiles = glob.glob(sdir+'/raw/*jp2')
        poolraw = Pool(processes=nproc)
        poolraw.map(remove_sym,rawfiles)
        poolraw.close()
        print 'Directories Already Exist. Proceeding to Download'
    
    
    #create a rough starting time for the weekly movie
    #assumes span is in days
    sday = eday+dt(days=-span)
    
    
    
    if addinfo:
        #HEK client to find number of events over the observation window
        hekc = hek.HEKClient()
        AR = 'AR'
        FL = 'FL'
        GOES = 'GOES'
        #add AR region labels
        arres = hekc.query(hek.attrs.Time(sday,eday),hek.attrs.EventType(AR))
        #find all GOES flares in time frame
        flres = hekc.query(hek.attrs.Time(sday,eday),hek.attrs.EventType(FL),hek.attrs.OBS.Instrument == GOES)
        
        #list of regions associated with the flare
        flreglist = []
        flintlist = []
        for i in flres:
            flreglist.append(i['ar_noaanum'])
        #just grab first letter for GOES class
            flintlist.append(i['fl_goescls'][0])
        #make flare list a numpy array so can run where  
        flreglist = np.array(flreglist)
        flintlist = np.array(flintlist)
        
        #unique active regions
        namelist = []
        #use to keep AR regions in order
        objelist = []
        #counting ARs
        for i in arres:
            arname = i['ar_noaanum']
        #check that the active region is named and is not previously added to the list (preventing duplicates)
            if ((arname is not None) & (arname not in namelist)):
                namelist.append(arname)
                objelist.append(i)
        
        # sort objects by AR number
        sorter   = np.argsort(namelist)
        namelist = np.array(namelist)[sorter]
        objelist = np.array(objelist)[sorter]
        
        #Classes of flares to count
        fltype = ['X','M','C']
        #create a dictionary containing the AR object, flare count, and 
        eventobj = {}
        for j,i in enumerate(objelist):
            ar = namelist[j] 
            flnum, = np.where(flreglist == ar) 
        #count the intensities of flares 
            arfla = flintlist[flnum]
            flarc = {}
            for gclass in fltype:
                foundg, = np.where(arfla == gclass)
                flarc[gclass] = foundg.size
        #store flare counts in dictionary
            eventobj[str(ar)+'_sflare'] = 'FL = {3:3d} (X={0:2d},M={1:2d},C={2:2d})'.format(flarc['X'],flarc['M'],flarc['C'],flnum.size) 
        #now find the end time of the event    
            endtimes = []
            for tar in arres:
                if tar['ar_noaanum'] == ar : endtimes.append(datetime.strptime(tar['event_endtime'],'%Y-%m-%dT%H:%M:%S')) 
            eventobj[str(ar)+'_endtime'] = max(endtimes)
#            print 'AR {0:5d} endtime = {1}'.format(ar,eventobj[str(ar)+'_endtime'].strftime('%Y-%m-%dT%H:%M:%S'))
    else:
        objelist = []
    
    #find times of previous files
    findp = glob.glob(sdir+'/working/*png')
    #array which will carry the cadence
    dayarray = []
    #if a new directory then begin with a three day time line
    if len(findp) == 0:
        sday = sday
        dayarray.append(sday.strftime('%Y/%m/%d %H:%M:%S'))
    #if files already exist build upon existing cadence
    else:
        for date in findp:
    # cut file into date array
            cdate = date.replace(sdir+'/working/','').replace('__SDO_AIA_AIA_193.png','')[:-3]
    # format into day
            mdate = datetime.strptime(cdate,"%Y_%m_%d__%H_%M_%S")
    #remove files greater than 3 days old
            if mdate < sday:
               print 'Removing {0} because it is old'.format(date)
               os.remove(date)
    #files are in numerical order so when it is done stop 
            else:
               ind = mdate.strftime('%Y/%m/%d %H:%M:%S') 
               dayarray.append(ind) 
    #build the start date from the start of the current list of files if files already exist 
    #find times of remaining previous files
        findp = glob.glob(sdir+'/working/*png')
    #convert file string into a day
    #This is useful when hv has a bad day and data is not pipelined correctly
    #build cadence from last available day
        date = findp[-1]
        cdate = date.replace(sdir+'/working/','').replace('__SDO_AIA_AIA_193.png','')[:-3]
        mdate = datetime.strptime(cdate,"%Y_%m_%d__%H_%M_%S")
        sday = mdate 
    
    #write ffmpeg command 
    com = open(sdir+'/run_ffmpeg.csh','w')
    com.write('/usr/local/bin/ffmpeg -y -f image2 -r {2:2d} -i working/symlinks/seq%4d.png -an -pix_fmt "yuv420p" -vcodec libx264 -level 41 -crf 18.0 -b 8192k -r {2:2d} -bufsize 8192k -maxrate 8192k -g 25 -coder 1 -profile main -preset faster -qdiff 4 -qcomp 0.7 -directpred 3 -flags +loop+mv4 -cmp +chroma -partitions +parti4x4+partp8x8+partb8x8 -subq 7 -me_range 16 -keyint_min 1 -sc_threshold 40 -i_qfactor 0.71 -rc_eq "blurCplx^(1-qComp)" -s "{0:4d}x{1:4d}" -b_strategy 1 -bidir_refine 1 -refs 6 -deblockalpha 0 -deblockbeta 0 -trellis 1 -x264opts keyint=25:min-keyint=1:bframes=1 -threads 2 working/{3}.mp4\n'.format(w0,h0,frate,outmov))
    com.close()
    
    
    
    
    #the first day is the start day
    nday = sday
    i = len(dayarray)
    j = 1
    samples = round(minweek/cadence*(span/7.))#number of dt samples to probe
#    while ((i < samples) & (nday < now-dt(hours=dh,minutes=dm))):
    while (nday < now-dt(hours=dh,minutes=dm)):
    #while ((i < 1) & (nday < now-dt(hours=dh))):
        nday = sday+dt(minutes=j*cadence) 
    #format the string for input
        ind = nday.strftime('%Y/%m/%d %H:%M:%S')
        dayarray.append(ind)
        i += 1
        j += 1 
    
    
    
    #for testing
    #dayarray = dayarray[:2]
    #download jp2 from hv
    pool = Pool(processes=nproc)
    outs = pool.map(get_file,dayarray)
    pool.close()
    
    #get a list of files
    dayarrayn = glob.glob(sdir+'/raw/*jp2')
    dayarraya = glob.glob(sdir+'/working/*png')
    for i in dayarrayn:
        dayarraya.append(i)
    
    #write new file and add text for time
    #create png files and symbolic links
    forpool = np.arange(len(dayarraya))
    pool1 = Pool(processes=nproc)
    outs = pool1.map(format_img,forpool)
    pool1.close()
    
    #change to current directory
    os.chdir(sdir)
    
    #change file to executable
    mod = subprocess.call(['/bin/csh','-c','chmod a+x run_ffmpeg.csh'])
    
    #run ffmpeg
    run = subprocess.call(['/bin/csh','-c','./run_ffmpeg.csh'])
    #rename the file to the file directory (cuts down on down time if running in a loop)
    os.rename('{1}/working/{0}.mp4'.format(outmov,sdir),'{1}/final/{0}.mp4'.format(outmov,sdir))


if __name__ == "__main__":
    main(sys.argv[1:])

