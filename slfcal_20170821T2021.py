from suncasa.utils import helioimage2fits as hf
import os
import numpy as np
import pickle
import pdb
import mod_slftbs as mods
from matplotlib import gridspec as gridspec
from sunpy import map as smap
from matplotlib import pyplot as plt
import time

#task handles
dofullsun=0
domasks=0
doslfcal=0
dofinalclean=0
doapply=0
doclean_slfcaled=1
workdir='/srg/bchen/EOVSA/solar/20170821_2020flare/'
slfcaldir=workdir+'slfcal_2021/'
imagedir=slfcaldir+'images/'
maskdir=slfcaldir+'masks/'
imagedir_final=slfcaldir+'images_final/'
#imagedir_slfcaled=slfcaldir+'images_slfcaled_155418/'
imagedir_slfcaled=slfcaldir+'images_slfcaled_2021/'
imagedir_slfcaled2=slfcaldir+'images_slfcaled_202120/'
caltbdir=slfcaldir+'caltbs/'
if not os.path.exists(imagedir):
    os.makedirs(imagedir)
if not os.path.exists(maskdir):
    os.makedirs(maskdir)
if not os.path.exists(imagedir_final):
    os.makedirs(imagedir_final)
if not os.path.exists(imagedir_slfcaled):
    os.makedirs(imagedir_slfcaled)
if not os.path.exists(imagedir_slfcaled2):
    os.makedirs(imagedir_slfcaled2)
if not os.path.exists(caltbdir):
    os.makedirs(caltbdir)
#refms = workdir+'IDB20170910163625.corrected.ms'
refms = workdir+'msdata/IDB20170821201020-203020.corrected.ms'
msapply = workdir+'msdata/IDB20170821201020-203020.12s.corrected.ms'
msapply_slfcaled = workdir+'msdata/IDB20170821201020-203020.12s.slfcaled.ms'
refms_slfcal = slfcaldir + os.path.basename(refms) + '.T202112-202114.XX.slfcal'
refms_slfcaled = slfcaldir + os.path.basename(refms) + '.T202112-202114.XX.slfcaled'
#refms_slfcal = refms + '.slfcal'
#refms_slfcaled = refms + '.slfcaled'
if not os.path.exists(refms_slfcal):
    #os.system('cp -r '+refms+' '+refms_slfcal)
    split(vis=refms,outputvis=refms_slfcal,timerange='20:21:12~20:21:14',datacolumn='data',correlation='XX')
# selected times for generating self-calibration solutions
tranges=['2017/08/21/20:21:12~2017/08/21/20:21:14']
spws=[str(s+1) for s in range(30)]
antennas='0~12' 
npix=512
nround=3 #number of slfcal cycles
# parameters specific to the event (found from imaging)
phasecenter='J2000 10h02m59 11d58m14'
xran=[280,480]
yran=[-50,150]

subms_a=[]
slfcalms_a=[]
slfcaledms_a=[]
for t,trange in enumerate(tranges):
    subms_ = slfcaldir+'IDB20171020.ms.t{0:d}.slfcal'.format(t)
    slfcalms_ = slfcaldir+'IDB20171020.ms.t{0:d}.XX.slfcal'.format(t)
    slfcaledms_ = slfcaldir+'IDB20171020.ms.t{0:d}.XX.slfcaled'.format(t)
    if not os.path.exists(slfcalms_):
        split(vis=refms,outputvis=slfcalms_,datacolumn='data',timerange=trange,correlation='XX')
    if not os.path.exists(subms_):
        split(vis=refms,outputvis=subms_,datacolumn='data',timerange=trange,correlation='')
    slfcalms_a.append(slfcalms_)
    subms_a.append(subms_)
    slfcaledms_a.append(slfcaledms_)

if dofullsun:
    #initial mfs clean to find out the image phase center
    im_init='tstimg_init'
    os.system('rm -rf '+im_init+'*')
    clean(vis=refms,
            antenna='0~12',
            imagename=im_init,
            spw='1~15',
            mode='mfs',
            timerange=trange,
            imagermode='csclean',
            psfmode='clark',
            imsize=[npix],
            cell=['5arcsec'],
            niter=1000,
            gain=0.05,
            stokes='I',
            restoringbeam=['30arcsec'],
            interactive=False,
            pbcor=True,
            usescratch=True)

    hf.imreg(vis=refms,imagefile=im_init+'.image',fitsfile=im_init+'.fits',
             timerange=trange,usephacenter=False,verbose=True)
    clnjunks = ['.flux', '.mask', '.model', '.psf', '.residual']
    for clnjunk in clnjunks:
        if os.path.exists(im_init + clnjunk):
            shutil.rmtree(im_init + clnjunk)

    from sunpy import map as smap
    from matplotlib import pyplot as plt
    eomap=smap.Map(im_init+'.fits')
    eomap.data=eomap.data.reshape((npix,npix))
    eomap.plot_settings['cmap'] = plt.get_cmap('jet')
    eomap.plot()
    eomap.draw_limb()
    eomap.draw_grid()
    plt.show()
    viewer(im_init+'.image')

masks_a=[]
if domasks:
    for t,trange in enumerate(tranges):
        slfcalms=slfcalms_a[t]
        clearcal(slfcalms)
        delmod(slfcalms)
        antennas='0~12' 
        pol='I'
        calprefix=caltbdir+'slf_t{0:d}'.format(t)
        imgprefix=maskdir+'slf_t{0:d}'.format(t)

        # step 1: set up the clean masks
        img_init=imgprefix+'_init_ar_'
        os.system('rm -rf '+img_init+'*')
        spwrans=['1~5','6~12','13~20','21~30']
        mask=[]
        imnames=[]
        for spwran in spwrans:
            imname=img_init+spwran.replace('~','-')
            try:
                clean(vis=slfcalms,
                        antenna='0~12',
                        imagename=imname,
                        spw=spwran,
                        mode='mfs',
                        timerange=trange,
                        imagermode='csclean',
                        psfmode='clark',
                        imsize=[npix],
                        cell=['2arcsec'],
                        niter=1000,
                        gain=0.05,
                        stokes='XX',
                        restoringbeam=['20arcsec'],
                        phasecenter=phasecenter,
                        weighting='briggs',
                        robust=1.0,
                        interactive=True,
                        pbcor=True,
                        usescratch=False)
                imnames.append(imname+'.image')
                mask.append(imname+'.mask')
                clnjunks = ['.flux', '.model', '.psf', '.residual']
                for clnjunk in clnjunks:
                    if os.path.exists(imname + clnjunk):
                        shutil.rmtree(imname + clnjunk)
            except:
                print('error in cleaning spw: '+spwran)
        masks_a.append(mask)

    pickle.dump(masks_a,open(slfcaldir+'masks_a.p','wb'))
    #viewer(imnames)

    #os.system('rm -rf '+calprefix+'*')
    #os.system('rm -rf '+imgprefix+'.spw*')
if doslfcal:
    if not os.path.exists(slfcaldir+'masks_a.p') or not ('masks_a' in vars()):
        print 'masks do not exist. Please run dopartsun first!'
        masks_a=[]
    if os.path.exists(slfcaldir+'masks_a.p'):
        masks_a=pickle.load(open(slfcaldir+'masks_a.p','rb'))
    os.system('rm -rf '+imagedir+'*')
    #first step: make a mock caltable for the entire database
    slftbs_a=[]
    for t,trange in enumerate(tranges):
        print('Processing '+str(t+1)+' in '+str(len(tranges))+' times: '+trange)
        slftbs=[]
        slfcalms=slfcalms_a[t]
        slfcaledms=slfcaledms_a[t]
        subms=subms_a[t]
        masks=masks_a[t]
        calprefix=caltbdir+'slf_t{0:d}'.format(t)
        imgprefix=imagedir+'slf_t{0:d}'.format(t)
        tb.open(slfcalms+'/SPECTRAL_WINDOW')
        reffreqs=tb.getcol('REF_FREQUENCY')
        bdwds=tb.getcol('TOTAL_BANDWIDTH')
        cfreqs=reffreqs+bdwds/2.
        tb.close()
        sbeam=40.
        strtmp=[m.replace(':','') for m in trange.split('~')]
        timestr='t'+strtmp[0]+'-'+strtmp[1]
        refantenna='0'
        niters=[100,300,500]
        robusts=[1.0,0.5,0.0]
        doapplycal=[1,1,1]
        calmodes=['p','p','a']
        uvranges=['','',''] 
        spwrans=['1~5','6~12','13~20','21~30']
        for n in range(nround):
            slfcal_tb_g= calprefix+'.G'+str(n)
            fig = plt.figure(figsize=(12,10))
            gs = gridspec.GridSpec(6, 5)
            for s,sp in enumerate(spws):
                print 'processing spw: '+sp
                cfreq=cfreqs[int(sp)]
                bm=max(sbeam*cfreqs[1]/cfreq,10.)
                slfcal_img = imgprefix+'.spw'+sp.zfill(2)+'.slfcal'+str(n)
                if n == 0:
                    spbg=max(int(sp)-2,1)
                    sped=min(int(sp)+2,30)
                    spwran=str(spbg)+'~'+str(sped)
                    if int(sp) <= 5:
                        mask = masks[0]
                    if int(sp) > 5 and int(sp) <= 12:
                        mask = masks[1]
                    if int(sp) > 12 and int(sp) <= 20:
                        mask = masks[2]
                    if int(sp) > 20 and int(sp) <= 30:
                        mask = masks[3]
                    print 'using spw {0:s} as model'.format(spwran)
                    print 'using mask {0:s}'.format(mask)
                else:
                    spwran = sp
                    mask = masks[0]
                try:
                    clean(vis=slfcalms,
                            antenna=antennas,
                            imagename=slfcal_img,
                            uvrange=uvranges[n],
                            #spw=sp,
                            spw=spwran,
                            mode='mfs',
                            timerange=trange,
                            imagermode='csclean',
                            psfmode='clark',
                            imsize=[npix],
                            cell=['2arcsec'],
                            niter=niters[n],
                            gain=0.05,
                            #stokes='I',
                            stokes='XX', #use pol XX image as the model
                            weighting='briggs',
                            robust=robusts[n],
                            phasecenter=phasecenter,
                            #mask='box [ [ 75pix , 90pix] , [205pix, 165pix ] ]',
                            mask=mask,
                            restoringbeam=[str(bm)+'arcsec'],
                            pbcor=False,
                            interactive=False,
                            usescratch=True)
                    if os.path.exists(slfcal_img+'.image'):
                        fitsfile=slfcal_img+'.fits'
                        hf.imreg(vis=slfcalms,imagefile=slfcal_img+'.image',fitsfile=fitsfile,
                                 timerange=trange,usephacenter=False,toTb=True,verbose=False,overwrite=True)
                    clnjunks = ['.mask','.flux', '.model', '.psf', '.residual', '.image']
                    for clnjunk in clnjunks:
                        if os.path.exists(slfcal_img + clnjunk):
                            shutil.rmtree(slfcal_img + clnjunk)
                    ax = fig.add_subplot(gs[s])
                    eomap=smap.Map(fitsfile)
                    eomap.data=eomap.data.reshape((npix,npix))
                    eomap.plot_settings['cmap'] = plt.get_cmap('jet')
                    eomap.plot()
                    eomap.draw_limb()
                    eomap.draw_grid()
                    ax.set_title(' ')
                    ax.set_xlim(xran)
                    ax.set_ylim(yran)

                except:
                    print 'error in cleaning spw: '+sp
                    print 'using nearby spws for initial model'
                    sp_e=int(sp)+2
                    sp_i=int(sp)-2
                    if sp_i < 0:
                        sp_i = 0
                    if sp_e > 30:
                        sp_e = 30
                    sp_=str(sp_i)+'~'+str(sp_e)
                    try:
                        tget(clean)
                        spw=sp_
                        clean()
                    except:
                        print 'still not successful. abort...'
                        break

                # copy model from xx to yy
                #pdb.set_trace()
                #mods.cpxx2yy(slfcalms)
                gaincal(vis=slfcalms, refant=refantenna,antenna=antennas,caltable=slfcal_tb_g,spw=sp, uvrange='',\
                        gaintable=[],selectdata=True,timerange=trange,solint='inf',gaintype='G',calmode=calmodes[n],\
                        combine='',minblperant=4,minsnr=2,append=True)
                if not os.path.exists(slfcal_tb_g):
                    #listcal(vis=slfcalms, caltable=slfcal_table)
                    print 'No solution found in spw: '+sp
            figname=imagedir+'slf_t{0:d}_n{1:d}.png'.format(t,n)
            plt.savefig(figname)
            time.sleep(10)
            plt.close()

            if os.path.exists(slfcal_tb_g):
                slftbs.append(slfcal_tb_g)
                slftb=[slfcal_tb_g]
                os.chdir(slfcaldir)
                if calmodes[n] == 'p': 
                    plotcal(caltable=slfcal_tb_g,antenna='1~12',xaxis='freq',yaxis='phase',\
                            subplot=431,plotrange=[-1,-1,-180,180],iteration='antenna',figfile=slfcal_tb_g+'.png',showgui=False)
                if calmodes[n] == 'a':
                    plotcal(caltable=slfcal_tb_g,antenna='1~12',xaxis='freq',yaxis='amp',\
                            subplot=431,plotrange=[-1,-1,0,2.],iteration='antenna',figfile=slfcal_tb_g+'.png',showgui=False)
                os.chdir(workdir)


            if doapplycal[n]:
                clearcal(slfcalms)
                delmod(slfcalms)
                applycal(vis=slfcalms,gaintable=slftb,spw=','.join(spws),selectdata=True,\
                         antenna=antennas,interp='nearest',flagbackup=False,applymode='calonly',calwt=False)

            if n < nround-1: 
                prompt=raw_input('Continuing to selfcal?')
                #prompt='y'
                if prompt.lower() == 'n':
                    if os.path.exists(slfcaledms):
                        os.system('rm -rf '+slfcaledms)
                    split(subms,slfcaledms,datacolumn='corrected')
                    print 'Final calibrated ms is {0:s}'.format(slfcaledms)
                    break
                if prompt.lower() == 'y':
                    slfcalms_=slfcalms+str(n)
                    if os.path.exists(slfcalms_):
                        os.system('rm -rf '+slfcalms_)
                    split(slfcalms,slfcalms_,datacolumn='corrected')
                    slfcalms=slfcalms_
            else:
                if os.path.exists(slfcaledms):
                    os.system('rm -rf '+slfcaledms)
                split(slfcalms,slfcaledms,datacolumn='corrected')
                print 'Final calibrated ms is {0:s}'.format(slfcaledms)

        slftbs_a.append(slftbs)

if dofinalclean:
    import glob
    pol='XX'
    if os.path.exists(slfcaldir+'masks_a.p'):
        masks_a=pickle.load(open(slfcaldir+'masks_a.p','rb'))
    for t,trange in enumerate(tranges):
        print('Processing '+str(t+1)+' in '+str(len(tranges))+' times: '+trange)
        img_final=imagedir_final+'/slf_final_{0}_t{1:d}'.format(pol,t)
        masks=masks_a[t]
        slfcaledms=slfcaledms_a[t]
        spws=[str(s+1) for s in range(30)]
        tb.open(slfcaledms+'/SPECTRAL_WINDOW')
        reffreqs=tb.getcol('REF_FREQUENCY')
        bdwds=tb.getcol('TOTAL_BANDWIDTH')
        cfreqs=reffreqs+bdwds/2.
        tb.close()
        sbeam=35.
        fig = plt.figure(figsize=(12,10))
        gs = gridspec.GridSpec(6, 5)
        for s,sp in enumerate(spws):
            cfreq=cfreqs[int(sp)]
            bm=max(sbeam*cfreqs[1]/cfreq,10.)
            imname=img_final+'_s'+sp.zfill(2)
            fitsfile=imname+'.fits'
            if int(sp) < 5:
                mask = masks[0]
            if int(sp) > 5 and int(sp) <= 12:
                mask = masks[1]
            if int(sp) > 12 and int(sp) <= 20:
                mask = masks[2]
            if int(sp) > 20 and int(sp) <= 30:
                mask = masks[3]
            #print 'using mask {0:s}'.format(mask)
            if not os.path.exists(fitsfile):
                print 'cleaning spw {0:s} with beam size {1:.1f}"'.format(sp,bm)
                try:
                    clean(vis=slfcaledms,
                            antenna=antennas,
                            imagename=imname,
                            spw=sp,
                            #mode='channel',
                            mode='mfs',
                            timerange=trange,
                            imagermode='csclean',
                            psfmode='clark',
                            imsize=[npix],
                            cell=['2arcsec'],
                            niter=1000,
                            gain=0.05,
                            stokes=pol,
                            weighting='briggs',
                            robust=0.0,
                            restoringbeam=[str(bm)+'arcsec'],
                            phasecenter=phasecenter,
                            #mask=mask,
                            mask='',
                            pbcor=True,
                            interactive=False,
                            usescratch=False)
                except:
                    print 'cleaning spw '+sp+' unsuccessful. Proceed to next spw'
                    continue
                if os.path.exists(imname+'.image'):
                    hf.imreg(vis=slfcaledms,imagefile=imname+'.image',fitsfile=fitsfile,
                             timerange=trange,usephacenter=False,toTb=True,verbose=False)
                junks=['.flux','.model','.psf','.residual','.mask','.image']
                for junk in junks:
                    if os.path.exists(imname+junk):
                        os.system('rm -rf '+imname+junk)
            else:
                ax = fig.add_subplot(gs[s])
                eomap=smap.Map(fitsfile)
                eomap.data=eomap.data.reshape((npix,npix))
                eomap.plot_settings['cmap'] = plt.get_cmap('jet')
                eomap.plot()
                eomap.draw_limb()
                eomap.draw_grid()
                ax.set_title(' ')
                ax.set_xlim(xran)
                ax.set_ylim(yran)
        figname=imagedir+'slf_t{0:d}_final.png'.format(t)
        plt.savefig(figname)
        plt.close()

if doapply:
    import glob
    slftbs_comb=[]
    for n in range(nround):
        tbin=glob.glob(caltbdir+'slf_t*.G{0}'.format(n))
        mods.concat(tb_in=tbin,tb_out=caltbdir+'slf_comb.G{0}'.format(n))
        slftbs_comb.append(caltbdir+'slf_comb.G{0}'.format(n))
    os.chdir(slfcaldir)
    plotcal(caltable=slftbs_comb[0], antenna='1~12',spw='1~10',xaxis='time',yaxis='phase',subplot=341,
            iteration='antenna',plotrange=[-1,-1,-100,100])
    os.chdir(workdir)
    clearcal(refms_slfcal)
    applycal(vis=refms_slfcal,gaintable=slftbs_comb,spw=','.join(spws),selectdata=True,\
             antenna=antennas,interp='linear',flagbackup=False,applymode='calonly',calwt=False)
    applycal(vis=msapply,gaintable=slftbs_comb,spw=','.join(spws),selectdata=True,\
             antenna=antennas,interp='linear',flagbackup=False,applymode='calonly',calwt=False)
    if os.path.exists(refms_slfcaled):
        os.system('rm -rf '+refms_slfcaled)
    split(refms_slfcal,refms_slfcaled,datacolumn='corrected')
    split(msapply,msapply_slfcaled,datacolumn='corrected')

if doclean_slfcaled:
    import glob
    pol='XX'
    if os.path.exists(slfcaldir+'masks_a.p'):
        masks_a=pickle.load(open(slfcaldir+'masks_a.p','rb'))
    #for t,trange in enumerate(tranges):
    for t in range(1):
        #trange=tranges[0]
        trange='2017/08/21/20:21:10~2017/08/21/20:21:30'
        print('Processing '+str(t+1)+' in '+str(len(tranges))+' times: '+trange)
        #img_final=imagedir_slfcaled+'/slf_final_{0}_t{1:d}'.format(pol,t)
        #vis = refms_slfcaled
        img_final=imagedir_slfcaled2+'/slf_final_{0}_t{1:d}'.format(pol,t)
        vis = msapply_slfcaled
        masks=masks_a[t]
        #spws=[str(s+1) for s in range(12)]
        spws=[str(s+1) for s in range(30)]
        tb.open(vis+'/SPECTRAL_WINDOW')
        reffreqs=tb.getcol('REF_FREQUENCY')
        bdwds=tb.getcol('TOTAL_BANDWIDTH')
        cfreqs=reffreqs+bdwds/2.
        tb.close()
        sbeam=30.
        from matplotlib import gridspec as gridspec
        from sunpy import map as smap
        from matplotlib import pyplot as plt
        fig = plt.figure(figsize=(12,10))
        gs = gridspec.GridSpec(5, 6)
        #fig = plt.figure(figsize=(12,9))
        #gs = gridspec.GridSpec(3, 4)
        for s,sp in enumerate(spws):
            cfreq=cfreqs[int(sp)]
            bm=max(sbeam*cfreqs[1]/cfreq,4.)
            imname=img_final+'_s'+sp.zfill(2)
            fitsfile=imname+'.fits'
            if s < 4:
                mask=''
            else:
                mask='circle [[281pix, 259pix] ,90pix ]'
            if not os.path.exists(fitsfile):
                print 'cleaning spw {0:s} with beam size {1:.1f}"'.format(sp,bm)
                try:
                    clean(vis=vis,
                            antenna=antennas,
                            imagename=imname,
                            spw=sp,
                            #mode='channel',
                            mode='mfs',
                            timerange=trange,
                            imagermode='csclean',
                            psfmode='clark',
                            imsize=[npix],
                            cell=['2arcsec'],
                            niter=1000,
                            gain=0.05,
                            stokes=pol,
                            weighting='briggs',
                            robust=2.0,
                            restoringbeam=[str(bm)+'arcsec'],
                            #restoringbeam=[],
                            phasecenter=phasecenter,
                            mask='',
                            pbcor=True,
                            interactive=False,
                            usescratch=False)
                except:
                    print 'cleaning spw '+sp+' unsuccessful. Proceed to next spw'
                    continue
                if os.path.exists(imname+'.image'):
                    hf.imreg(vis=vis,imagefile=imname+'.image',fitsfile=fitsfile,
                             timerange=trange,usephacenter=False,toTb=True,verbose=False)
                junks=['.flux','.model','.psf','.residual','.mask','.image']
                for junk in junks:
                    if os.path.exists(imname+junk):
                        os.system('rm -rf '+imname+junk)
            else:
                ax = fig.add_subplot(gs[s])
                eomap=smap.Map(fitsfile)
                eomap.data=eomap.data.reshape((npix,npix))
                eomap.plot_settings['cmap'] = plt.get_cmap('jet')
                eomap.plot()
                eomap.draw_limb()
                eomap.draw_grid()
                ax.set_title(' ')
                ax.get_xaxis().set_visible(False)
                ax.get_yaxis().set_visible(False)
                ax.set_xlim(xran)
                ax.set_ylim(yran)
                plt.text(0.98,0.9,'{0:.1f} GHz'.format(cfreq/1e9),transform=ax.transAxes,ha='right',color='w',fontweight='bold')
        plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
        figname=imagedir_slfcaled2+'slf_202120_final.png'
        plt.savefig(figname)
        plt.close()

