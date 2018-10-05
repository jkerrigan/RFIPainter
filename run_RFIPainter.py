import matplotlib
matplotlib.use('TKAgg')
import pylab as plt
from scipy import signal
import aipy as a
import numpy as np
from glob import glob
import random
import pyuvdata
from matplotlib.patches import Rectangle
import h5py
import RFIPainter_tools as rpt

"""
RFI Painter allows you to import miriad or hdf5 files containing waterfall visibilities
to flag RFI or genenerally any bad data.
Important keys to note are:

Left Click: place an RFI flag
Left Click + Hold and Drag: place a rectanglar grouping of RFI flags
Right Click: undo the most recently placed set of flags
p: do not include waterfall into output dataset
right arrow: go to next waterfall
q: quit out of the canvasses, do this when you've reached the end of your dataset
   and all of the waterfalls that were 'right arrowed' will be saved
o: enter zoom mode

"""

pull_samples=False
#file_locs = '/Users/josh/Desktop/Organize/UPennML/ml_rfi/RealVisRFI_v3.h5'
file_locs = '/Users/josh/Desktop/RFIMLDiverseDataset/'
suffix = 'uvSLIM'
output_name = 'RawRealVis_v1.h5'

plt.ion()

mask = None
# Load data either from miriad files or hdf5 files
# Miriad file requires a suffix (e.g. uv or uvOCRS)

data,mask,info = rpt.loadPYUVdata(file_locs,suffix)
#data,mask = rpt.loadH5(file_locs)

info = {}
print('Input data size: {} (# Samples, times, freqs)'.format(np.shape(data)))
paint_session = rpt.Painter(data,info,mask=mask)
plt.show(block=True)
data,mask,updated_info = paint_session.get_data()

data = [data[int(i.split('_')[-1])] for i in updated_info.keys() if int(updated_info[i].split('_')[0]) != -1]
mask = [mask[int(i.split('_')[-1])] for i in updated_info.keys() if int(updated_info[i].split('_')[0]) != -1]
print('Shape of data: {0} Shape of mask: {1}'.format(np.shape(data),np.shape(mask)))

try:
    updated_info_str = ['{}_'.format(i.split('/')[1].split('_')[0])+updated_info[i] for i in updated_info.keys() if int(updated_info[i].split('_')[0]) != -1]
except:
    updated_info_str = []
print(updated_info_str)

mask_ = np.copy(mask) 
mask = np.logical_not(mask)

hf = h5py.File(output_name,'a') # caution! this will attempt to overwrite data
suffix='uv'
if pull_samples:
    for suffix in suffixes:
        try:
            dp = hf.create_group(suffix)
        except:
            dp = hf.require_group(suffix)
        HERAdata = []
        new_mask = []
        obs_info = []
        for key in updated_info.keys():
            print key
            obs = ('.').join(key.split('.')[:-1])+'.'+suffix+'SLIM'
            uv = pyuvdata.UVData()
            uv.read_miriad(obs)
            antpairs = np.array(uv.get_antpairs())
            order = updated_info[key]
            rnd_aps = antpairs[random.sample(range(0,len(antpairs)),99)]
            for ap in rnd_aps:
                data_grab = uv.get_data(ap).squeeze()
                refresher = Refresh(data_grab,info,mask=mask_[order],delay_only=True)
                plt.show(block=True)
                keep = refresher.return_keep()
                if keep: 
                    if len(HERAdata) == 0:
                        HERAdata = [data_grab]
                        new_mask = [mask[order]]
                        obs_info = [obs+'_'+str(ap)]
                    else:
                        HERAdata.append(data_grab)
                        new_mask.append(mask[order])
                        obs_info.append(obs+'_'+str(ap))
                else:
                    print('Waterfall rejected.')
                    continue
        sh = np.shape(HERAdata)
        try:
            dp.create_dataset('data', data=HERAdata)
            dp.create_dataset('flag', data=new_mask)
            dp.create_dataset('meta', data=obs_info)
        except:
            old_data = np.copy(dp['data'][:])
            old_flags = np.copy(dp['flag'][:])
            old_info = np.copy(dp['meta'][:])
            hf.close()
            hf = h5py.File(output_name, 'w')
            dp = hf.create_group(suffix)
            new_data = np.vstack((old_data,HERAdata))
            new_flags = np.vstack((old_flags,new_mask))
            new_info = np.vstack((old_info,obs_info))
            print(np.shape(new_data))
            dp.create_dataset('data', data=new_data)
            dp.create_dataset('flag', data=new_flags)
            dp.create_dataset('meta', data=new_info)
    hf.close()
else:
    try:
        updated_info_str = ['{}_'.format(i.split('/')[1])+updated_info[i] for i in updated_info.keys()]
        print(updated_info_str)
    except:
        pass
    try:
        dp = hf.create_group(suffix)
    except:
        dp = hf.require_group(suffix)
    try:
        dp.create_dataset('data', data=data)
        dp.create_dataset('flag', data=mask)
        dp.create_dataset('meta', data=updated_info_str)
    except:
        old_data = np.copy(dp['data'][:])
        old_flags = np.copy(dp['flag'][:])
        old_info = np.copy(dp['meta'][:])
        hf.close()
        hf = h5py.File(output_name, 'w')
        dp = hf.create_group(suffix)
        new_data = np.vstack((old_data,data))
        new_flags = np.vstack((old_flags,mask))
        new_info = np.append(old_info,updated_info_str)
        print(np.shape(new_data))
        dp.create_dataset('data', data=new_data)
        dp.create_dataset('flag', data=new_flags)
        dp.create_dataset('meta', data=new_info)
    hf.close()

