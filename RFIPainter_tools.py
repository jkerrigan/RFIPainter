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

Plot=False
pull_samples=False

def loadH5(file_locs):
    f = h5py.File(file_locs)
    try:
        return f['data'][:],np.logical_not(f['flag'])
    except:
        return f['uv']['data'][:],np.logical_not(f['uv']['flag'])
def loadPYUVdata(file_locs,suffix):
    files = np.array(glob(file_locs+'*'+suffix))
    info = {}
    ct = 0
    rnd_ = np.random.randint(len(files),size=5)
    file_cut = files[rnd_]
    for f in file_cut:
        print(f)
        uv = pyuvdata.UVData()
        uv.read_miriad(f,run_check_acceptability=False,run_check=False,check_extra=False)
        antpairs = uv.get_antpairs()
        flag_npz_name = '.'.join(f.split('.')[:5])+'.uvOC.flags.npz'
        try:
            flag_npz = np.load(flag_npz_name)
        except:
            #continue
            print('Skipping because npz flag file -{}- not found'.format(flag_npz_name))
        for i in range(10):
            rnd = np.random.randint(len(antpairs))
            ap1,ap2 = antpairs[rnd]
            bsl = uv.antnums_to_baseline(ap1,ap2)
            try:
                HERAlabels_ = np.logical_not(flag_npz['flag_array'][bsl == flag_npz['baseline_array']])
            except:
                HERAlabels_ = 0.
                #pass
            pos1 = uv.antenna_positions[uv.antenna_numbers==ap1]
            pos2 = uv.antenna_positions[uv.antenna_numbers==ap2]
            bl_len = np.round(np.sqrt(np.sum((pos1-pos2)**2)),2)
            info[f+'_{0}'.format(ct)] = '{0}_blen_{1}'.format(ct,bl_len)#str(ct)+'_blen_'+str() #['antpairs'][ct] = antpairs[rnd] #antpairs[rnd] : ct}
            if f == file_cut[0] and i == 0:
                HERAdata = [uv.get_data(antpairs[rnd]).squeeze()]
                try:
                    HERAlabels = [HERAlabels_.squeeze()]
                except:
                    HERAlabels = []
            else:
                HERAdata.append(uv.get_data(antpairs[rnd]).squeeze())
                try:
                    HERAlabels.append(HERAlabels_.squeeze())
                except:
                    HERAlabels = []
            ct+=1
        del(uv)
    print('Dataset size: ',np.shape(HERAdata))
    if np.dim(HERAlabels) > 1:
        HERAlabels = np.zeros_like(HERAdata).real
    return HERAdata,HERAlabels,info

plt.ion()
options={'y':-1,
    'n':100}

class Annotate(object):
    def __init__(self):
        self.ax = plt.gca()
#        self.ax = 
        self.rect = Rectangle((0,0), 1, 1)
        self.x0 = None
        self.y0 = None
        self.x1 = None
        self.y1 = None
        self.ax.add_patch(self.rect)
        self.ax.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.ax.figure.canvas.mpl_connect('button_release_event', self.on_release)
        self.ax.figure.canvas.mpl_connect('key_press_event', self.quit)

    def on_press(self, event):
        #print 'press'
        self.x0 = event.xdata
        self.y0 = event.ydata

    def on_release(self, event):
        #print 'release'
        self.x1 = event.xdata
        self.y1 = event.ydata
        self.rect.set_width(self.x1 - self.x0)
        self.rect.set_height(self.y1 - self.y0)
        self.rect.set_xy((self.x0, self.y0))
        self.ax.figure.canvas.draw()

class Painter:
    def __init__(self,data,info,mask=None,delay_only=False):
        self.ct = 0
        if delay_only:
            self.fig2 = plt.figure()
            self.canvas2 = self.fig2.canvas
        else:
            self.fig = plt.figure()
            self.fig2 = plt.figure()
            self.fig3 = plt.figure()
            self.canvas = self.fig.canvas
            self.canvas2 = self.fig2.canvas
            self.canvas3 = self.fig3.canvas
        self.data = np.copy(data)
        self.keep = True
        self.delay_only = delay_only
        if np.any(mask!=None):
            self.mask = mask
        else:
            self.mask = np.ones_like(self.data).astype(int)
        self.info = info
        self.curr_mask = self.mask[0]
        self.undomask = np.copy(self.curr_mask)
        self.curr_data = self.data[0]
        if delay_only:
            self.ax2 = self.fig2.gca()
        else:
            self.ax = self.fig.gca()
            self.ax2 = self.fig2.gca()
            self.ax3 = self.fig3.gca()
        if delay_only:
            self.cid7_2 = self.canvas2.mpl_connect('key_press_event', self.pop)
        else:
            self.cid1 = self.canvas.mpl_connect('button_press_event', self.on_press)
            self.cid2 = self.canvas.mpl_connect('button_release_event', self.on_release)
            self.cid3 = self.canvas.mpl_connect('button_press_event', self.undo)
            self.cid4 = self.canvas.mpl_connect('key_press_event', self.on_key_press)
            self.cid5 = self.canvas.mpl_connect('key_release_event', self.on_key_release)
            self.cid6 = self.canvas.mpl_connect('key_press_event', self.next)
            self.cid7 = self.canvas.mpl_connect('key_press_event', self.pop)
        
            self.cid1_ = self.canvas3.mpl_connect('button_press_event', self.on_press)
            self.cid2_ = self.canvas3.mpl_connect('button_release_event', self.on_release)
            self.cid3_ = self.canvas3.mpl_connect('button_press_event', self.undo)
            self.cid4_ = self.canvas3.mpl_connect('key_press_event', self.on_key_press)
            self.cid5_ = self.canvas3.mpl_connect('key_release_event', self.on_key_release)
            self.cid6_ = self.canvas3.mpl_connect('key_press_event', self.next)
            self.cid7_ = self.canvas3.mpl_connect('key_press_event', self.pop)

        logdata = np.log10(np.abs(self.curr_data))*np.abs(self.curr_mask)

        maxD = np.max(np.abs(logdata))
        minD = np.min(np.abs(logdata))
        if self.delay_only:
            self.delay_plot = self.ax2.imshow(np.log10(np.abs(self.delayTrans())),aspect='auto')
            self.ax2.set_title('Delay Spectrum')
            self.ax2.set_xlabel('Delay Bin')
            self.ax2.set_ylabel('Time')
        else:
            self.amp_plot = self.ax.imshow(np.log10(np.abs(self.curr_data))*np.abs(self.curr_mask),aspect='auto')
            self.ax.set_title('Visibility')
            self.ax.set_xlabel('Freq.')
            self.ax.set_ylabel('Time')
            self.delay_plot = self.ax2.imshow(np.log10(np.abs(self.delayTrans())),aspect='auto')
            self.ax2.set_title('Delay Spectrum')
            self.ax2.set_xlabel('Delay Bin')
            self.ax2.set_ylabel('Time')
            self.phs_plot = self.ax3.imshow(np.angle(self.curr_data*self.curr_mask),aspect='auto')
            self.ax3.set_title('Phase')
            self.ax3.set_xlabel('Freq')
            self.ax3.set_ylabel('Time')
        self.x0 = None
        self.y0 = None
        self.x1 = None
        self.y1 = None

    def on_key_press(self, event):
        if event.key == 'shift':
            self.shift_is_held = True

    def test(self, event):
        if event.key == 'right':
            print('Right Arrow')
        if event.key == 'n':
            print('n')

    def pop(self, event):
        if event.key == 'p':
            if self.delay_only:
                self.keep = False
            else:
                if self.info.has_key(self.ct):
                    bad_key = [key_ for key_ in self.info.keys() if int(info[key_].split('_')[0]) == self.ct]
                    self.info[bad_key[0]] = '{0}_blen_{1}'.format(-1,000)
                else:
                    self.info['000_'+str(self.ct)] = '{0}_blen_{1}'.format(-1,000)
                self.ct += 1
                self.curr_data = self.data[self.ct]
                self.curr_mask = self.mask[self.ct]
                if self.delay_only:
                    self.delay_plot.set_data(np.log10(np.abs(self.delayTrans())))
                    self.canvas2.draw()
                else:
                    self.amp_plot.set_data(np.log10(np.abs(self.curr_data))*np.abs(self.curr_mask))
                    self.canvas.draw()
                    self.delay_plot.set_data(np.log10(np.abs(self.delayTrans())))
                    self.canvas2.draw()
                    self.phs_plot.set_data(np.angle(self.curr_data*self.curr_mask))
                    self.canvas3.draw()
    
    def next(self, event):
        if event.key == 'right':
            self.data[self.ct] = np.copy(self.curr_data)
            self.mask[self.ct] = np.copy(self.curr_mask)
            if not self.info.has_key(self.ct):
                self.info['000_'+str(self.ct)] = '{0}_blen_{1}'.format(self.ct,000)
            self.ct += 1
            print(self.ct)
            self.curr_data = self.data[self.ct]
            self.curr_mask = self.mask[self.ct]
            if self.delay_only:
                self.delay_plot.set_data(np.log10(np.abs(self.delayTrans())))
                self.canvas2.draw()
            else:
                self.amp_plot.set_data(np.log10(np.abs(self.curr_data))*np.abs(self.curr_mask))
                self.canvas.draw()
                self.delay_plot.set_data(np.log10(np.abs(self.delayTrans())))
                self.canvas2.draw()
                self.phs_plot.set_data(np.angle(self.curr_data*self.curr_mask))
                self.canvas3.draw()

    def on_key_release(self, event):
        if event.key == 'shift':
            self.shift_is_held = False

    def onclick(self,event):
        print('button=%d, x=%d, y=%d, xdata=%f, ydata=%f' % (event.button, event.x, event.y, event.xdata, event.ydata))
        self.undomask = self.mask
        self.curr_mask[int(event.ydata)+1,int(event.xdata)+1] = 0

    def undo(self, event):
        if event.button == 3:
            self.curr_mask = self.undomask
            logdata = np.log10(np.abs(self.curr_data))*np.abs(self.curr_mask)
            minD = np.min(logdata)
            maxD = np.max(logdata)
            if self.delay_only:
                self.delay_plot.set_data(np.log10(np.abs(self.delayTrans())))
                self.canvas2.draw()
            else:
                self.amp_plot.set_data(np.log10(np.abs(self.curr_data))*np.abs(self.curr_mask))
                self.canvas.draw()
                self.delay_plot.set_data(np.log10(np.abs(self.delayTrans())))
                self.canvas2.draw()
                self.phs_plot.set_data(np.angle(self.curr_data*self.curr_mask))
                self.canvas3.draw()

    def on_press(self, event):
        if event.button == 1:
            self.undomask = np.copy(self.curr_mask)
            self.x0 = event.xdata
            self.y0 = event.ydata

    def quit(self, event):
        if event.button == 'q':
            self.canvas.mpl_disconnect(self.cid1)
            self.canvas.mpl_disconnect(self.cid2)
            self.canvas.mpl_disconnect(self.cid3)
            self.canvas.mpl_disconnect(self.cid4)
            self.canvas.mpl_disconnect(self.cid5)
            self.ax.close()
            self.ax2.close()
            

    def on_release(self, event):
        if plt.get_current_fig_manager().toolbar.mode != '': return

        if event.button == 1:
            self.x1 = event.xdata
            self.y1 = event.ydata
            if self.y0 >= self.y1 and self.x0 >= self.x1:
                self.curr_mask[int(self.y1):int(self.y0+1),int(self.x1):int(self.x0+1)] = 0
            elif self.y0 >= self.y1 and self.x1 >= self.x0:
                self.curr_mask[int(self.y1):int(self.y0+1),int(self.x0):int(self.x1+1)] = 0
            elif self.x0 >= self.x1 and self.y1 >= self.y0:
                self.curr_mask[int(self.y0):int(self.y1+1),int(self.x0):int(self.x1+1)] = 0
            else:
                self.curr_mask[int(self.y0):int(self.y1+1),int(self.x0):int(self.x1+1)] = 0
            logdata = np.log10(np.abs(self.curr_data))*np.abs(self.curr_mask)
            minD = np.min(logdata)
            maxD = np.max(logdata)
            if self.delay_only:
                self.delay_plot.set_data(np.log10(np.abs(self.delayTrans())))
                self.canvas2.draw()
            else:
                self.amp_plot.set_data(np.log10(np.abs(self.curr_data))*np.abs(self.curr_mask))
                self.canvas.draw()
                self.delay_plot.set_data(np.log10(np.abs(self.delayTrans())))
                self.canvas2.draw()
        
                self.phs_plot.set_data(np.angle(self.curr_data*self.curr_mask))
                self.canvas3.draw()

    
    def get_data(self):
        return self.data,self.mask,self.info

    def return_keep(self):
        return self.keep

    def delayTrans(self):
        bh = a.dsp.gen_window(1024,window='blackman-harris')
        DATA = np.fft.fft(self.curr_data*np.abs(self.curr_mask)*bh,axis=1)
        DATA_ = np.fft.fftshift(DATA,axes=1)
        return DATA_
