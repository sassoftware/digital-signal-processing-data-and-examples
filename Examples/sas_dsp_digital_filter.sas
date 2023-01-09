/*********************************************************************************************************
Digital Filter Example

This example illustrates how to design a digital filter and apply it to 
time series data using digital filtering functions in SAS/IML.

For more information about digital signal processing functions in SAS/IML, 
please refer to: 
https://go.documentation.sas.com/doc/en/pgmsascdc/v_034/imlug/imlug_dsp_toc.htm

Copyright © 2023, SAS Institute Inc., Cary, NC, USA.  All Rights Reserved.
SPDX-License-Identifier: Apache-2.0

*********************************************************************************************************/

/* 1. Generate the input signal */
proc iml;
	sf  = 1000;   /* sampling frequency */
	sf2 = sf/2;   /* Nyquist frequency = half the sampling frequency */
	T = 1.0/sf;   /* sampling period */
	N = 200;
	t = do(0, (N-1)*T, T);

	f1 = 5;
	f2 = 100;
	f3 = 200;
	pi = constant("pi");
	x1 = sin(2*pi*f1*t);
	x2 = sin(2*pi*f2*t);
	x3 = sin(2*pi*f3*t);
	x = 0.5*x1 + x2 + x3; /* input signal */

	create input_sig var {"t", "x"};
	append;
	close input_sig;
quit

title "Input of the Digital Filter";  
proc sgplot data=input_sig;
	series x=t y=x;
    xaxis grid label="Time";
    yaxis grid label="Amplitude";
run;

/* 2. Filter design and filter review */
proc iml;
	/* Fiter design */
	filter_name = "butter";
	filter_type = "bandpass";

	Wp = 0.15 || 0.25;
	Ws = 0.05 || 0.35;
	Rp = 3;
	Rs = 50;
	call dforder(n, Wc, filter_name, filter_type, Wp, Ws, Rp, Rs);
	print n, Wc;

	call dfdesign(b, a, z, p, k, filter_name, filter_type, n, Wc);
	print b a, z, p, k;  /* b a are printed together */

	/*Filter review*/
	N = 1024;
	h = DFSOSFreqzZPK(z, p, k, N);

	w = h[,3] / constant('pi');
	logh = 20*log10(sqrt(h[,1]#h[,1]+h[,2]#h[,2]) + 1e-16);
	/* a small constant is added to avoid the input to log10() is 0 */

	uc_x1 = do(-1, 1, 0.001);
	uc_x2 = do(1, -1, -0.001);
	uc_y1 = sqrt(1-uc_x1#uc_x1);
	uc_y2 = -sqrt(1-uc_x2#uc_x2);
	uc_x = uc_x1 || uc_x2;
	uc_y = uc_y1 || uc_y2;

	zero_re = z[, 1];
	zero_im = z[, 2];
	pole_re = p[, 1];
	pole_im = p[, 2];

	create filter_params var {"zero_re", "zero_im",
	                          "pole_re", "pole_im", "k"};
	append;
	close filter_params;

	create filter_plot var {"w", "logh", "zero_re", "zero_im", "pole_re",
		                    "pole_im", "uc_x", "uc_y"};
	append;
	close filter_plot;
quit;

/* Plot the filter's frequency response */
proc sgplot data=filter_plot aspect=0.667;
   title "Frequency Response";
   series x=w y=logh;
   xaxis grid label="Normalized frequency ((*ESC*){unicode pi} rad/sample)";
   yaxis grid label="Magnitude (dB)";
run;


/* 3. Apply the digital filter to input signal */
proc iml;

	use input_sig;
	read all var {"t" "x"};
	close input_sig;

	use filter_params;
	read all var {"zero_re" "zero_im" "pole_re" "pole_im" "k"};
	close filter_params;

	z = zero_re || zero_im;
	p = pole_re || pole_im;
	k = k[1];

	y = dfsosfilt(x, z, p, k);

	create filter_results var {"t", "x", "y"};
	append;
	close filter_results;
quit;

/* Plot the filterd signal */
proc sgplot data=filter_results;
	title "Input and Output of the Digital Filter";
	series x=t y=x / legendlabel="Input signal";
	series x=t y=y / legendlabel="Output signal" lineattrs=(thickness=2);
	xaxis grid label="Time";
	yaxis grid label="Amplitude";
run;
