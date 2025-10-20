import os 
import sys
from os import path

import pickle
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from itertools import product
from time import time 
from datetime import datetime
#import seaborn as sns
from glob import glob

import ROOT
import uproot

from keras.models import Sequential, Model
from keras.layers import Dense, Input, Dropout, BatchNormalization
from keras.utils import plot_model
from keras.callbacks import EarlyStopping, Callback, ReduceLROnPlateau, ModelCheckpoint
from keras.constraints import unit_norm
from keras.optimizers import SGD, Adam

import sklearn as sk
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve
from sklearn.preprocessing import RobustScaler, StandardScaler

sys.path.append('../')
#sys.path.append('../../objects')
from tools import Tools
from mva_tools import MVATools
#from samples import signal_samples
#from baseline_selection import selection
#from categories import categories

#TODO try with CNN?
#TODO use phi_k1_dztrg etc. as feature

def getOptions():
  from argparse import ArgumentParser
  parser = ArgumentParser(description='Script to perform the NN training', add_help=True)
  parser.add_argument('--category_batch', type=str, dest='category_batch', help='category label'                               , default=None)
  parser.add_argument('--outdir'        , type=str, dest='outdir'        , help='output directory'                             , default=None)
  parser.add_argument('--process'       , dest='process'                 , help='run the process function', action='store_true', default=False)
  return parser.parse_args()


class Sample(object):
  '''
    Class to convert the sample into dataframe while applying some selection
  '''
  def __init__(self, filename, selection):
    self.filename = filename
    self.selection = selection
    self.mva_tools = MVATools()
    try:
        with uproot.open(self.filename) as file:
          tree = file['signal_tree']
          df = tree.arrays(library="pd")
          self.df = df.query(self.mva_tools.getPandasQuery(self.selection))
    except:
        print('No entry was found with the requested selection')
        self.df = pd.DataFrame()
    print('\t selected events: {}'.format(len(self.df)))


class Trainer(object):
  #def __init__(self, features, epochs, batch_size, learning_rate, number_nodes, scaler_type, do_early_stopping, do_reduce_lr, signal_label, data_pl, data_tagnano, data_tagflat, nsigma, dirname, baseline_selection, categories, category_batch, outdir):
  def __init__(self, features, epochs, batch_size, learning_rate, number_nodes, scaler_type, do_early_stopping, do_reduce_lr, signal_label, data_pl, data_tagnano, data_tagflat, nsigma, dirname, baseline_selection, outdir):
    self.mva_tools = MVATools()

    self.features = features
    self.epochs = epochs
    self.batch_size = batch_size
    self.learning_rate = learning_rate
    self.number_nodes = number_nodes
    self.scaler_type = scaler_type
    self.do_early_stopping = do_early_stopping
    self.do_reduce_lr = do_reduce_lr
    self.dirname = dirname + '_' + datetime.now().strftime('%Y%b%d_%Hh%Mm%Ss')

    self.baseline_selection = baseline_selection
    #self.categories = categories
    #self.category_batch = category_batch
    self.signal_label = signal_label
    #self.signal_files = signal_samples[self.signal_label]

    self.outdir = outdir

    self.username = 'anlyon'
    self.data_pl = data_pl
    self.data_tagnano = data_tagnano
    self.data_tagflat = data_tagflat

    self.nsigma = nsigma

    self.target_branch = 'is_signal'
    self.do_scale_key = True

    self.resolution_p0 = 6.98338e-04
    self.resolution_p1 = 7.78382e-03


  def createOutDir(self):
    '''
      This function creates the output directory
    '''
    outdir = './outputs/{}'.format(self.dirname)
    if not path.exists(outdir):
      os.system('mkdir -p {}'.format(outdir))

    return outdir


  def writeSubmitter(self, category):
    '''
      Write bash submitter
    '''
    content = '\n'.join([
        '#!/bin/bash',
        'homedir="$PWD"',
        'workdir=/scratch/{}/training_{}'.format(self.username, category.label),
        'mkdir -p $workdir',
        'cp -r ../*py $workdir', #TODO adapt if changing directory
        'cp -r ./*py $workdir',
        'cp -r ../../objects/*py $workdir',
        'cd $workdir',
        'DATE_START=`date +%s`',
        'python trainer.py --category_batch {cat} --outdir {out} --process'.format(cat=category.label, out=self.outdir),
        'DATE_END=`date +%s`',
        'runtime=$((DATE_END-DATE_START))',
        'echo " --> Wallclock running time: $runtime s"',
        'cp -r *h5 $homedir/{out}'.format(out=self.outdir),
        'cp -r *pck $homedir/{out}'.format(out=self.outdir),
        'cp -r *png $homedir/{out}'.format(out=self.outdir),
        'cp -r *pdf $homedir/{out}'.format(out=self.outdir),
        'cp -r *txt $homedir/{out}'.format(out=self.outdir),
        'cp trainer.py $homedir/{out}'.format(out=self.outdir),
        'cd $homedir',
        'rm -r $workdir',
        ])
    submitter_name = '{}/submitter_{}.sh'.format(self.outdir, category.label)
    submitter = open(submitter_name, 'w+')
    submitter.write(content)
    submitter.close()
    '\n -> {} created'.format(submitter_name)


  def submit(self, category):
    '''
      Submit bash script on slurm
    '''
    command = 'sbatch -p standard --account t3 -o {out}/log_{cat}.txt -e {out}/log_{cat}.txt --job-name=trainer_{cat} --mem {mem} {out}/submitter_{cat}.sh'.format(
        out = self.outdir,
        cat = category.label,
        mem = 20000,
        )
    print('\n --> submitting category {}'.format(category.label))
    os.system(command)


  def saveFig(self, plt, name):
    '''
      Save python figure
    '''
    plt.savefig('{}/{}.pdf'.format(self.outdir, name))    
    plt.savefig('{}/{}.png'.format(self.outdir, name))    
    print(' --> {}/{}.png created'.format(self.outdir, name))


  def getDataSamples(self, extra_selection=None, max_files=-1):
    '''
      Function that fetches the samples into lists
    '''
    #print('========> starting reading the trees'))
    now = time()

    ### data 
    #path = '/pnfs/psi.ch/cms/trivcat/store/user/{}/BHNLsGen/data'.format(self.username)
    #if self.data_tagflat != None:
    #  filename = '{path}/{pl}/ParkingBPH1_Run2018D/Chunk*/flat/flat_bparknano_{tagnano}_{tagflat}.root'.format(path=path, pl=self.data_pl, tagnano=self.data_tagnano, tagflat=self.data_tagflat)
    #else:
    #  filename = '{path}/{pl}/ParkingBPH1_Run2018D/Chunk*/flat/flat_bparknano_{tagnano}.root'.format(path=path, pl=self.data_pl, tagnano=self.data_tagnano)

    #FIXME hardcoded
    #glob_path = '/eos/user/a/anlyon/CPVGen/V00_D1/*/*/*/*0/flat'
    #glob_path = '/eos/cms/store/group/phys_bphys/anlyon/CPVGen/data/V03_24Sep25/2018/D1/ParkingBPH1/crab_V03_24Sep25_D1_20251002_201533/251002_181540/0000/flat'
    glob_path = '/eos/cms/store/group/phys_bphys/anlyon/CPVGen/data/V03_24Sep25/2018/A1/ParkingBPH1/crab_V03_24Sep25_A1_20250930_181836/250930_161914/0001/flat'
    filename = 'flat_bparknano.root'
    data_filenames = [f for f in glob('{}/{}'.format(glob_path, filename))]

    #data_filenames = sorted(data_filenames, key=lambda file_: float(file_[file_.find('Chunk')+5:file_.find('_', file_.find('Chunk')+5)]))

    data_samples = []
    for ifile, data_filename in enumerate(data_filenames):
      if max_files != -1 and ifile > max_files: continue
      print(data_filename)
      data_samples.append(Sample(filename=data_filename, selection=self.baseline_selection + ' && ' + extra_selection))

    print('========> it took %.2f seconds' %(time() - now))

    return data_samples


  def getMCSamples(self, extra_selection=None, max_files=-1):
    '''
      Function that fetches the samples into lists
    '''
    print('========> starting reading the trees')
    now = time()
    #FIXME hardcoded
    #filename_mc_1 = '/eos/user/a/anlyon/CPVGen/102X_crab_trgmu_filter/BsToPhiPhiTo4K/nanoFiles/merged/flat_bparknano.root'
    filename_mc_1 = '/eos/cms/store/group/phys_bphys/anlyon/CPVGen/signal_central/V03_24Sep25/2018/dgneq0/BsToPhiPhiTo4K_TuneCP5_13TeV_pythia8-evtgen/crab_V03_24Sep25_dgneq0_20251001_220102/251001_200142/0000/flat/flat_bparknano.root' 
    mc_samples = [
      Sample(filename=filename_mc_1, selection=self.baseline_selection + ' && ' + extra_selection),
    ]
      
    print('========> it took %.2f seconds' %(time() - now))

    return mc_samples


  def getMassList(self, is_bc=False):
    '''
      Get the list of signal masses used in the training
    '''
    masses = []
    for signal_file in self.signal_files:
      mass = signal_file.mass
      if mass not in masses:
        if not is_bc and mass <= 4.5:
          masses.append(mass)
        elif is_bc and mass >= 3.:
          masses.append(mass)

    masses.sort()

    return masses


  def removeInfs(self, df):
    '''
      Remove rows with infs/nan in at least one of the features
    '''
    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.dropna(subset=self.features)

    return df


  def createDataframe(self, data_samples, mc_samples):
    '''
      Function that converts root samples to pandas dataframe
    '''
    data_df = pd.concat([idt.df for idt in data_samples], sort=False)
    mc_df   = pd.concat([imc.df for imc in mc_samples]  , sort=False)

    data_df = self.removeInfs(data_df)
    mc_df = self.removeInfs(mc_df)
        
    data_df['weight'] = 1. #FIXME keep here?
    mc_df['weight'] = 1. #FIXME keep here? reweight mc here for balanced training?

    return data_df, mc_df


  def createParametrisedDataframe(self, extra_selection, data_type, is_bc=False, statistics=None, max_files=-1):
    '''
      Function to create the dataframe with the training features and parameter
      The statistics is balanced between the signal and background, among the signal
      hypotheses and the background windows
    '''
    if data_type not in ['data', 'mc']:
      raise RuntimeError('Unknown data type "{}". Please check'.format(data_type))

    if data_type == 'data' and statistics == None:
      raise RuntimeError('Please first create mc dataframe in order to get the available statistics')

    pd.options.mode.chained_assignment = None

    masses = self.getMassList(is_bc=is_bc)
    
    if data_type == 'mc':
      stats = dict.fromkeys(masses, 0)
      n_ctaus = dict.fromkeys(masses, 0)
      dfs = dict.fromkeys(masses, []) 

      for signal_file in self.signal_files: 
        # get sample
        signal_filename = signal_file.filename if not is_bc else signal_file.filename_Bc
        if is_bc and signal_file.filename_Bc == None: continue # skip points for which there is no Bc sample
        if not is_bc and signal_file.filename == None: continue # skip mass 5.5 GeV for non-Bc category
        sample = Sample(filename=signal_filename, selection=self.baseline_selection + ' && ' + extra_selection)

        # get dataframe
        the_df = sample.df

        # get the statistics per mass point
        stats[signal_file.mass] += len(the_df)

        # get number of ctau samples per mass point
        n_ctaus[signal_file.mass] += 1

        # set the mass parameter to the exact value
        the_df['mass_key'] = signal_file.mass

        # weight the mc with the trigger scale factors
        #the_df = the_df.rename(columns={'weight_hlt_D1': 'weight'})
        the_df['weight'] = 1.

        dfs[signal_file.mass] = dfs[signal_file.mass] + [the_df]

      # get the minimum statistics 
      #statistics = min(stats[min(stats, key=stats.get)], 5000)
      #aimed_statistics = stats[min(stats, key=stats.get)]
      statistics = stats[min(stats, key=stats.get)]
      
      df = pd.DataFrame()
      for mass in masses:
        df_tmp = pd.concat([(idt.sample(statistics/n_ctaus[mass]) if statistics/n_ctaus[mass]<len(idt) else idt) for idt in dfs[mass]], sort=False)
        df = pd.concat([df, df_tmp], sort=False) 
      
        
    elif data_type == 'data':
      # get the sample
      samples = self.getDataSamples(extra_selection=extra_selection, max_files=max_files) 

      df_full = [] # will be used for plotting purposes later on
      dfs = dict.fromkeys(masses, []) 

      for sample in samples:
        # get the dataframe
        the_df = sample.df
        df_full.append(the_df)

        # first set an invalid parameter for the full spectrum. Ranges with this parameter will be removed from the training
        the_df['mass_key'] = -1 

        # and define the windows of nsigma around the mass hypothesis
        # set the mass parameter for a nsigma window range
        window_check = -1
        for mass in masses:
          resolution = self.resolution_p0 + self.resolution_p1 * mass
          window_min = mass - self.nsigma * resolution
          window_max = mass + self.nsigma * resolution
          # prevent windows to overlap
          #if window_check != -1 and window_min <= window_check:
          #  raise RuntimeError('There is an overlap of window around mass {} GeV with the previous window. Please adapt the window size.'.format(mass))

          window = 'hnl_mass > {} && hnl_mass < {}'.format(window_min, window_max)
          df_window = the_df.query(self.mva_tools.getPandasQuery(window))
          df_window['mass_key'] = mass

          # set the weight
          df_window['weight'] = 1.

          dfs[mass] = dfs[mass] + [df_window]

          window_check = window_max

      df_full = pd.concat([idt for idt in df_full], sort=False)

      df = pd.DataFrame()
      for mass in masses:
        requested_stat = statistics
        stat_permass = 0
        for the_df in dfs[mass]:
          stat_permass += len(the_df)
        if stat_permass < statistics:
          print('Requesting {} events from a sample of {} events. Please provide a larger input data sample'.format(statistics, stat_permass))
          requested_stat = stat_permass
        df_tmp = pd.concat([idt for idt in dfs[mass]], sort=False)
        df_tmp = df_tmp.sample(requested_stat)
        df = pd.concat([df, df_tmp], sort=False) 
    
    df = self.removeInfs(df)

    if data_type == 'mc':
      return df, statistics
    else:
      return df, df_full


  def assignTarget(self, df, branch, target):
    '''
      Add the target branch to the dataframes
    '''
    pd.options.mode.chained_assignment = None
    df[branch] = target

    return df


  def doScaling(self, X):
      '''
        Normalise the input features with a keras scaler 
      '''
      if self.scaler_type == 'robust':
        qt = RobustScaler()
      elif self.scaler_type == 'standard':
        qt = StandardScaler()
      else:
        raise RuntimeError('Unknown scaler "{}" - Aborting...'.format(self.scaler_type))

      qt.fit(X[self.features])
      xx = qt.transform(X[self.features])

      return xx, qt


  #def preprocessing(self, data_df, mc_df, label):
  def preprocessing(self, data_df, mc_df):
    '''
      Preprocessing of data before training/testing the NN
      This includes:
        - building the main_df
        - building the scaler
        - get the scaled features xx
        - get the target Y
    '''

    # concatenate the events and shuffle
    main_df = pd.concat([data_df, mc_df], sort=False)

    # re-index
    main_df.index = np.array(range(len(main_df)))

    # shuffle
    main_df = main_df.sample(frac=1, replace=False, random_state=1986) # of course, keep R's seed ;)

    # X and Y
    X = pd.DataFrame(main_df, columns=list(set(self.features)))
    Y = pd.DataFrame(main_df, columns=[self.target_branch])

    # scale the features
    # this is an important step!
    xx, qt = self.doScaling(X)

    # and save the scaler, which will have to be used throughout the full process, even at evaluation time
    #scaler_filename = '/'.join([self.outdir, 'input_tranformation_weighted_{}.pck'.format(label)])
    scaler_filename = '/'.join([self.outdir, 'input_tranformation_weighted.pck'])
    pickle.dump(qt,open(scaler_filename, 'wb'))
    print(' --> {} created'.format(scaler_filename))

    # save the exact list of features
    #features_filename = '/'.join([self.outdir, 'input_features_{}.pck'.format(label)])
    features_filename = '/'.join([self.outdir, 'input_features.pck'])
    pickle.dump(self.features, open(features_filename, 'wb' ))

    print(' --> {} created'.format(features_filename))

    return main_df, qt, xx, Y


  def defineModel(self):
    '''
      Define the NN
    '''
    #NOTE for the moment, everything is hardcoded

    activation = 'relu'
    
    # define the net
    n_input = len(self.features) 
    input  = Input((n_input,))
    layer  = Dense(self.number_nodes, activation=activation, name='dense1', kernel_constraint=unit_norm())(input)
    output = Dense(1 , activation='sigmoid' , name='output', )(layer)

    # Define outputs of your model
    model = Model(input, output)
    optimiser = Adam(learning_rate=self.learning_rate)
    model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['mae', 'acc'])
    
    print(model.summary())

    return model


  def defineCallbacks(self, patience_es, patience_lr):
    '''
      Define the callbacks
    '''
    # early stopping
    monitor = 'val_loss'
    es = EarlyStopping(monitor=monitor, mode='auto', verbose=1, patience=patience_es)
    
    # reduce learning rate when at plateau, fine search the minimum
    reduce_lr = ReduceLROnPlateau(monitor=monitor, mode='auto', factor=0.2, patience=patience_lr, min_lr=0.00001, cooldown=10, verbose=True)
    
    # save the model every now and then
    # kept only during excecution time and removed afterwards
    filepath = '/'.join([self.outdir, 'saved-model-{epoch:04d}_val_loss_{val_loss:.4f}_val_acc_{val_acc:.4f}.h5'])
    save_model = ModelCheckpoint(filepath, monitor='val_acc', verbose=1, save_best_only=True, save_weights_only=False, mode='auto')
    
    callbacks = [save_model]
    if self.do_early_stopping:
      callbacks.append(es)

    if self.do_reduce_lr:
      callbacks.append(reduce_lr)

    return callbacks


  def prepareInputs(self, xx, Y):
    '''
      Separate the main dataframe into training and validation sets
      Note: the input xx should arlready be scaled
    '''

    #NOTE this function is not adapted to training with weight; use prepareScaledInputs instead

    x_train, x_val, y_train, y_val = train_test_split(xx, Y, test_size=0.2, shuffle=True)
    
    # name the columns (selected in doScaling) consistently
    x_train = pd.DataFrame(x_train, columns=list(set(self.features))) # alternative to X_train[self.features[:]]
    x_val = pd.DataFrame(x_val, columns=list(set(self.features)))

    x_train = x_train.reset_index(drop=True)
    x_val = x_val.reset_index(drop=True)
    y_train = y_train.reset_index(drop=True)
    y_val = y_val.reset_index(drop=True)

    return x_train, x_val, y_train, y_val


  def prepareScaledInputs(self, main_df, Y, qt):
    '''
      Separate the main dataframe into training and validation sets
      The features are normalised according to the scaler qt

      NB: function not adapted to the parametric case
    '''

    x_train_tot, x_val_tot, y_train, y_val = train_test_split(main_df, Y, test_size=0.2, shuffle=True)

    # select training features
    x_train = x_train_tot[self.features]
    x_val = x_val_tot[self.features]

    # scale the features
    x_train = qt.transform(x_train)
    x_val = qt.transform(x_val)

    # select weight
    weight_train = x_train_tot['weight']
    weight_val = x_val_tot['weight']
    weight_train = weight_train.reset_index(drop=True)
    weight_val = weight_val.reset_index(drop=True)

    return x_train, x_val, y_train, y_val, weight_train, weight_val


  def train(self, model, x_train, y_train, x_val, y_val, weight_train, weight_val, callbacks):
    '''
      Perform the training
    '''
    history = model.fit(x_train, y_train, validation_data=(x_val, y_val, weight_val), sample_weight=weight_train, epochs=self.epochs, callbacks=callbacks, batch_size=self.batch_size, verbose=True)

    return history
    

  def predictScore(self, model, df):
    '''
      Return score with scaled input features
    '''
    x = pd.DataFrame(df, columns=self.features)

    # apply the scaler
    scaler_filename = '/'.join([self.outdir, 'input_tranformation_weighted.pck'])
    qt = pickle.load(open(scaler_filename, 'rb'))

    # if not scaling key
    xx = qt.transform(x[self.features])

    # predict
    score = model.predict(xx)

    return score


  #def saveModel(self, model, label):
  def saveModel(self, model):
    '''
      Save the model
    '''
    model_filename = '{}/net_model_weighted.h5'.format(self.outdir)
    model.save(model_filename)
    print(' --> {} created'.format(model_filename))


  def plotParametrisedMass(self, df, df_full=None, data_type='', label=''):
    '''
    '''
    if data_type not in ['data', 'mc']:
      raise RuntimeError('Unknown data type "{}". Please check'.format(data_type))

    masses = self.getMassList()
    bin_min = 0
    bin_max = 6.5
    nbins = 80.

    fig = plt.figure()
    ax = fig.add_subplot(111)
    #if not df_full.empty:
    if isinstance(df_full, pd.DataFrame):
      ax.hist(df_full['hnl_mass'], bins=np.arange(bin_min, bin_max, (bin_max-bin_min)/nbins), color='whitesmoke', label='unused')

    for mass in masses:
      resolution = self.resolution_p0 + self.resolution_p1 * mass
      window_min = mass - self.nsigma * resolution
      window_max = mass + self.nsigma * resolution
      window = 'hnl_mass > {} && hnl_mass < {}'.format(window_min, window_max)
      df_window = df.query(self.mva_tools.getPandasQuery(window))
      ax.hist(df_window['hnl_mass'], bins=np.arange(bin_min, bin_max, (bin_max-bin_min)/nbins), label='mass key = {}'.format(mass))

    ax.legend(loc='upper right',prop={'size': 12})
    ax.set_title(label, fontsize=20)
    ax.set_xlabel('Mass spectrum',fontsize=18)
    self.saveFig(fig, 'mass_{}_{}'.format(data_type, label))
    plt.clf()


  def plotLoss(self, history):
    '''
      Plot the loss for training and validation sets
    '''
    loss_train = history.history['loss']
    loss_val = history.history['val_loss']
    epochs_range = range(1, len(loss_train)+1)
    epochs = epochs_range
    plt.plot(epochs, loss_train, 'g', label='Training loss')
    plt.plot(epochs, loss_val, 'b', label='Validation loss')
    plt.title('Training and Validation Loss')
    plt.xlabel('Epochs')
    plt.ylabel('Loss')
    plt.legend()
    self.saveFig(plt, 'loss')
    plt.clf()


  def plotAccuracy(self, history):
    '''
      Plot the accuracy for training and validation sets
    '''
    acc_train = history.history['acc']
    acc_val = history.history['val_acc']
    epochs_range = range(1, len(acc_train)+1)
    epochs = epochs_range
    plt.plot(epochs, acc_train, 'g', label='Training accuracy')
    plt.plot(epochs, acc_val, 'b', label='Validation accuracy')
    plt.title('Training and Validation Accuracy')
    plt.xlabel('Epochs')
    plt.ylabel('Accuracy')
    plt.legend(loc='lower left')
    self.saveFig(plt, 'accuracy')
    plt.clf()


  def plotScore(self, model, mc_test_df, data_test_df):
    '''
      Plot the score distributions for signal and background
    '''
    # get the score for the test dataframe
    sig_score = self.predictScore(model, mc_test_df)
    bkg_score = self.predictScore(model, data_test_df)
    # add the score to the dataframes
    mc_test_df['score'] = sig_score
    data_test_df['score'] = bkg_score

    # plot the score distributions
    fig = plt.figure()
    ax = fig.add_subplot(111)
    bkg_score = [data_test_df['score']]
    bkg_name=['data-driven background']
    ax.hist(bkg_score, bins=np.arange(0,1.025,0.025), color='blue', alpha=0.8, label=bkg_name)
    ax.hist(mc_test_df['score'], bins=np.arange(0,1.025,0.025), color='darkorange', alpha=1, label='signal', histtype='step', linewidth=2)
    ax.legend(loc='upper left',prop={'size': 12})
    ax.set_title("Score distribution for testing set", fontsize=20)
    ax.set_xlabel('Score',fontsize=18)
    self.saveFig(fig, 'score')
    plt.clf()


  def plotROC(self, model, x_train, y_train, x_val, y_val):
    '''
      Plot the ROC curve
      Note that the x inputs are already scaled
    '''
    score_train = model.predict(x_train)
    fpr, tpr, wps = roc_curve(y_train, score_train) 

    plt.plot(fpr, tpr, label='train ROC')
    #print("AUC train",sk.metrics.auc(fpr,tpr))
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')

    score_val = model.predict(x_val)
    fpr, tpr, wps = roc_curve(y_val, score_val) 
    plt.plot(fpr, tpr, label='test ROC')
    plt.title('ROC')
    #print("AUC test",sk.metrics.auc(fpr,tpr))

    xy = [i*j for i,j in product([10.**i for i in range(-8, 0)], [1,2,4,8])]+[1]
    plt.plot(xy, xy, color='grey', linestyle='--')
    plt.yscale('linear')
    plt.legend()
    self.saveFig(plt, 'ROC')
    plt.clf()


  #def plotCorrelations(self, model, df, data_type, label):
  #  '''
  #    Plot the correlation matrix based on the training set
  #  '''
  #  if data_type not in ['data', 'mc']:
  #    raise RuntimeError('Unknown data_type "{}". Aborting'.format(data_type))

  #  # get the score for the test dataframe
  #  score = self.predictScore(model, df, label)

  #  # add the score to the dataframes
  #  df['score'] = score

  #  corr = df[self.features + ['score']].corr()
  #  
  #  # Set up the matplotlib figure
  #  f, ax = plt.subplots(figsize=(11, 9))
  #  
  #  # Generate a custom diverging colormap
  #  cmap = sns.diverging_palette(220, 10, as_cmap=True)

  #  # Draw the heatmap with the mask and correct aspect ratio
  #  g = sns.heatmap(corr, cmap=cmap, vmax=1., vmin=-1, center=0, annot=True, fmt='.2f',
  #                  square=True, linewidths=.8, cbar_kws={"shrink": .8})

  #  # rotate axis labels
  #  g.set_xticklabels(self.features+['score'], rotation='vertical')
  #  g.set_yticklabels(self.features+['score'], rotation='horizontal')

  #  plt.title('Linear Correlation Matrix - {}'.format(data_type))
  #  plt.tight_layout()
  #  self.saveFig(plt, 'correlations_{}_{}'.format(data_type, label))
  #  plt.clf()


  def plotKSTest(self, model, x_train, x_val, y_train, y_val, data_type):
    '''
      Plot the outcome of the Kolmogorov test
      Statistical test of compatibility in shape between two histograms
      Used to test the overfitting
    '''
    if data_type not in ['data', 'mc']:
      raise RuntimeError('Unknown data_type "{}". Aborting'.format(data_type))
    
    x_train = pd.DataFrame(x_train, columns=list(set(self.features)))
    x_val = pd.DataFrame(x_val, columns=list(set(self.features)))

    x_train = x_train.reset_index(drop=True)
    x_val = x_val.reset_index(drop=True)
    y_train = y_train.reset_index(drop=True)
    y_val = y_val.reset_index(drop=True)

    # only keep the data or mc components of the features
    # does it mess up with the normalisation?
    if data_type == 'data':
      x_train_part = x_train.drop(y_train.query('is_signal==1').index)
      x_val_part = x_val.drop(y_val.query('is_signal==1').index)
    else:
      x_train_part = x_train.drop(y_train.query('is_signal==0').index)
      x_val_part = x_val.drop(y_val.query('is_signal==0').index)

    score_train = model.predict(x_train_part)
    score_val = model.predict(x_val_part)

    h1 = ROOT.TH1F('train', 'train', 30, -1, 1)
    h2 = ROOT.TH1F('val', 'val', 30, -1, 1)
    for t, b in zip(score_train, score_val):
      h1.Fill(t)
      h2.Fill(b)

    c1=ROOT.TCanvas()
    if h1.Integral()!=0: h1.Scale(1./h1.Integral())
    if h2.Integral()!=0: h2.Scale(1./h2.Integral())
    c1.Draw()
    #h1.SetTitle(label)
    h1.Draw("hist")
    h2.SetLineColor(ROOT.kRed)
    h2.SetFillColor(ROOT.kWhite)
    h1.SetFillColor(ROOT.kWhite)
    h2.Draw("hist SAME")
    c1.BuildLegend()
    ks_score = h1.KolmogorovTest(h2)
    ks_value = ROOT.TPaveText(0.7, 0.65, 0.88, 0.72, 'nbNDC')
    ks_value.AddText('KS score {} = {}'.format(data_type, round(ks_score, 3)))
    ks_value.SetFillColor(0)
    ks_value.Draw('EP same')

    c1.SaveAs('{}/KS_test_{}.png'.format(self.outdir, data_type))


  def process(self):
    print('---- MVA Trainer ----')
    
    # create output directory
    if self.outdir == None:
      print('\n -> create output directory')
      self.outdir = self.createOutDir()
    else:
      self.outdir = '.'

    # copy script
    os.system('cp trainer.py {}'.format(self.outdir))
        
    #TODO use year as category
    #for category in self.categories:
    #if category.label == 'incl': continue
    #if category.label != 'lxysiggt150_OS': continue
    #if category.label != 'lxysig50to150_OS': continue
    #print('\n-.-.-')
    #print('category: {}'.format(category.label))
    #print('-.-.-')

    #if self.category_batch != None and category.label != category_batch: continue 

    # get the samples
    print('\n -> get the samples')
    data_samples = self.getDataSamples(extra_selection = 'abs(bs_mass_corr - 5.367) > 0.2')
    mc_samples = self.getMCSamples(extra_selection='ismatched == 1')

    # create dataframes
    print('\n -> create the dataframes')
    data_df, mc_df = self.createDataframe(data_samples, mc_samples)

    #else:
    #  print('\n -> create the dataframes')
    #  # do not load too much files in case of large statistics
    #  #if category.label in ['lxysig0to50_OS', 'lxysig0to50_SS', 'lxysig0to50_OS_Bc', 'lxysig0to50_SS_Bc']: max_files = 10 
    #  #elif category.label in ['lxysig50to150_OS', 'lxysig50to150_SS', 'lxysig50to150_OS_Bc', 'lxysig50to150_SS_Bc']: max_files = 30 
    #  #else: max_files = -1
    #  max_files = -1

    #  mc_df, statistics = self.createParametrisedDataframe(extra_selection=category.definition_flat, data_type='mc', is_bc=category.is_bc)
    #  data_df, data_df_full = self.createParametrisedDataframe(extra_selection=category.definition_flat, data_type='data', statistics=statistics, max_files=max_files, is_bc=category.is_bc)

    # assign the signal tag
    data_df = self.assignTarget(data_df, self.target_branch, 0) 
    mc_df = self.assignTarget(mc_df, self.target_branch, 1) 

    # preprocessing the dataframes
    print('\n -> preprocessing the dataframes' )
    #main_df, qt, xx, Y = self.preprocessing(data_df, mc_df, category.label)
    main_df, qt, xx, Y = self.preprocessing(data_df, mc_df)

    # define the NN
    print('\n -> defining the model' )
    model = self.defineModel()

    # define the callbacks
    print('\n -> defining the callbacks' )
    #if statistics < 3000:
    #  patience_es = 5
    #  patience_lr = 3
    #else:
    patience_es = 10
    patience_lr = 5
    #callbacks = self.defineCallbacks(category.label, patience_es, patience_lr)
    callbacks = self.defineCallbacks(patience_es, patience_lr)

    # out of the main_df, define which data chunks to 
    # train and test on. Make sure that the features
    # are scaled
    print('\n -> prepare the inputs' )
    x_train, x_val, y_train, y_val, weight_train, weight_val = self.prepareScaledInputs(main_df, Y, qt)
    #x_train, x_val, y_train, y_val = self.prepareInputs(xx, Y)

    ## create statistics file
    #stat_filename = '{}/statistics_{}.txt'.format(self.outdir, category.label)
    #stat_file = open(stat_filename, 'w+')
    #stat_file.write('\nAimed statistics per mass parameter: {}'.format(statistics))
    #stat_file.write('\nFull statistics for data: {}'.format(len(data_df)))
    #stat_file.write('\nFull statistics for mc: {}'.format(len(mc_df)))
    #stat_file.write('\nFull statistics for data+mc: {}'.format(len(main_df)))
    #stat_file.write('\nData+mc statistics used for training: {}'.format(len(x_train)))
    #stat_file.write('\nData+mc statistics used for validation: {}'.format(len(x_val)))
    #stat_file.close()
    #print(' --> {} created'.format(stat_filename) )

    # do the training
    print('\n -> training...' )
    history = self.train(model, x_train, y_train, x_val, y_val, weight_train, weight_val, callbacks)

    # save the model
    print('\n -> save the model' )
    #self.saveModel(model, category.label)
    self.saveModel(model)

    # plotting
    print('\n -> plotting...' )
    self.plotLoss(history)
    self.plotAccuracy(history)
    #self.plotScore(model, mc_test_df, data_test_df, category.label)
    #self.plotScore(model, mc_df, data_df, category.label)
    self.plotROC(model, x_train, y_train, x_val, y_val)
    #self.plotCorrelations(model, data_df, 'data', category.label)
    #self.plotCorrelations(model, mc_df, 'mc', category.label)
    self.plotKSTest(model, x_train, x_val, y_train, y_val, 'data')
    self.plotKSTest(model, x_train, x_val, y_train, y_val, 'mc')

    # cleaning
    print('\n -> cleaning')
    os.system('rm -r {}/saved-model*h5'.format(self.outdir))

    print('\n --- Done ---')


  def process_batch(self):
    print('---- MVA Trainer ----')

    # create output directory
    print('\n -> create output directory')
    self.outdir = self.createOutDir()

    for category in self.categories:
      if category.label == 'incl': continue
      if category.label != 'lxysig50to150_OS': continue
      print('\n-.-.-')
      print('category: {}'.format(category.label))
      print('-.-.-')
      self.writeSubmitter(category=category)
      self.submit(category=category)
      
    print('\n --- Done ---')



if __name__ == '__main__':
  ROOT.gROOT.SetBatch(True)

  #TODO try to balance the training, and to run on more data
  #https://en.innovatiana.com/post/how-to-balance-training-datasets

  features = [
    'bs_pt',
    'bs_eta',
    #'bs_beta',
    'bs_lxysig', 
    'bs_cos2d',
    'bs_sv_prob',
    'bs_mass_err',
    #'cos_theta_star_phi1',
    #'deltar_min',
    #'k1k3_mass',
    #'k1k4_mass',
    #'k2k3_mass',
    #'k2k4_mass',
    'phi1_pt', 
    'phi2_pt', 
    #'phi1_mass', 
    #'phi2_mass', 
    'phi1_mass_err', 
    'phi2_mass_err', 
    'phi1_cos2d',
    'phi2_cos2d',
    'phi1_sv_prob',
    'phi2_sv_prob',
    'phi1_lxysig',
    'phi2_lxysig',
    'k1_pt',
    'k2_pt',
    'k3_pt',
    'k4_pt',
    'k1_eta',
    'k2_eta',
    'k3_eta',
    'k4_eta',
    'k1_dcasig',
    'k2_dcasig',
    'k3_dcasig',
    'k4_dcasig',
    'k1_dzsig',
    'k2_dzsig',
    'k3_dzsig',
    'k4_dzsig',
    'k1_numberpixellayers',
    'k2_numberpixellayers',
    'k3_numberpixellayers',
    'k4_numberpixellayers',
    'k1_numbertrackerlayers',
    'k2_numbertrackerlayers',
    'k3_numbertrackerlayers',
    'k4_numbertrackerlayers',
    'k1_numberofvalidhits',
    'k2_numberofvalidhits',
    'k3_numberofvalidhits',
    'k4_numberofvalidhits',
  ]

  #features = ['pi_pt','mu_pt', 'mu0_pt','hnl_cos2d', 'sv_lxysig', 'pi_dcasig_corr', 'sv_prob', 'mu0_mu_mass', 'mu0_pi_mass', 'b_mass','deltar_mu0_mu', 'deltar_mu0_pi', 'mu0_pfiso03_rel', 'mu_pfiso03_rel', 'pi_numberofpixellayers', 'pi_numberoftrackerlayers', 'mu_numberofpixellayers', 'mu_numberoftrackerlayers', 'mu0_numberofpixellayers', 'mu0_numberoftrackerlayers']

  epochs = 3#50#3 #100
  batch_size = 50 #32 #50 #32
  learning_rate = 0.005 #0.01
  number_nodes = 64
  scaler_type = 'robust'
  do_early_stopping = True
  do_reduce_lr = True

  dirname = 'test'

  baseline_selection = 'bs_charge==0'

  #categories = categories['categories_0_50_150_Bc']
  #category_batch = getOptions().category_batch

  outdir = getOptions().outdir

  submit_batch = False

  nsigma = 10
  signal_label = 'V13_06Feb23_m2'
  data_pl = 'V13_06Feb23'
  data_tagnano = '06Feb23'
  data_tagflat = '31Jul23'

  trainer = Trainer(
      features = features, 
      epochs = epochs,
      batch_size = batch_size,
      learning_rate = learning_rate,
      number_nodes = number_nodes,
      scaler_type = scaler_type,
      do_early_stopping = do_early_stopping,
      do_reduce_lr = do_reduce_lr,
      signal_label = signal_label, 
      data_pl = data_pl,
      data_tagnano = data_tagnano,
      data_tagflat = data_tagflat,
      nsigma = nsigma,
      dirname = dirname,
      baseline_selection = baseline_selection,
      #categories = categories,
      #category_batch = category_batch,
      outdir = outdir,
      )

  if submit_batch and not getOptions().process:
    trainer.process_batch()
  else:
    trainer.process()




