/*********************************************************************************************************
Time Frequency Analysis Example - Siren Sound

In this example, the short-time Fourier transform is used to show how the 
frequency spectrum of an emergency siren changes over time. Autodetection 
of sirens at traffic intersections can alert drivers to make way for 
emergency vehicles and reduce collisions related to emergency vehicles. 
The sound of a siren shows a unique pattern in the frequency domain when
you use time-frequency analysis.

For more information about digital signal processing functions in SAS/IML, 
please refer to: 
https://go.documentation.sas.com/doc/en/pgmsascdc/v_034/imlug/imlug_dsp_toc.htm

Copyright © 2023, SAS Institute Inc., Cary, NC, USA.  All Rights Reserved.
SPDX-License-Identifier: Apache-2.0

*********************************************************************************************************/

/* Load the audio data 
   The sound data is saved in "siren.csv", which can be downloaded from:
   https://github.com/sassoftware/digital-signal-processing-data-and-examples/blob/main/Data/Audio_Files/siren.csv
*/ 

PROC IMPORT OUT= WORK.siren
DATAFILE= "siren.csv"
DBMS=CSV REPLACE;
GETNAMES=NO;
DATAROW=1;
RUN;

ods graphics /NXYBINSMAX= 2000000 MAXOBS=3000000;

proc iml;
   use work.siren;
   read all var _NUM_ into x;
   close work.siren;

   Fs  = 22050;         /* sampling frequency of the siren sound data*/
   win_len = 0.2 * Fs;  /* 200 ms in time duration */
   dt  = 1/Fs;
   n = nrow(x);
   t = do(1, n*dt, dt); /* timestamps of the data vector */

   overlap = round((win_len+1)/2);
   fftlen = 2##(ceil(log2(win_len)));
   window = tfwindow(win_len, "Hanning");

   /* 1. Use TFSTFT and SPECTROVIZ subroutines*/
   y = tfstft(x, window, overlap, fftlen);

   freq = y[, 1]*Fs;
   time = y[, 2]/Fs;
   dbPower = 10*log10(abs(y[, 3])+1e-12);

   labls = {"Time","Frequency","Power (dB)","Avg dBPower"};
   call spectroviz(t, x, dbPower, time, freq)
        panel=0
        labels=labls
        title="STFT Power Using TFSTFT and SPECTROVIZ"
        colorramp="BlueGreenRed";

   /* 2. Use SPECTROGRAM subroutine*/
   call spectrogram(x, win_len)
        fftlen=fftlen
        sampRate=Fs
        panel=0
        title="STFT Power Using SPECTROGRAM"
        colorramp="BlueGreenRed"
        scale='DB';

   /* Thresholding on the time-frequency heat map */
   thresholdDB = 35;
   T = dbPower > thresholdDB; /* Binarize the power values */
   dbPower_new = dbPower # T;
   labls = {"Time","Frequency","Power (dB)"};
   call spectroviz(t, x, dbPower_new, time, freq)
        panel=0
        labels=labls
        title="STFT Power after Thresholding"
        colorramp="BlueGreenRed";
quit;
