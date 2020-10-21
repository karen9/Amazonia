__Commands used in MrBayes__

```
prset brlenspr = clock:fossilization;
prset clockvarpr = igr;
prset igrvarpr = exp(10);
prset clockratepr=normal(0.01,0.005);

prset sampleprob = 0.7;
prset samplestrat = diversity;
prset speciationpr = exp(10);
prset extinctionpr = beta(1,1); 
prset fossilizationpr = beta(1,1); 
	
outgroup Callicebus;
constraint root = 1-.;
constraint rootdos  = 1-3;
constraint roottrees = 5-15;       
constraint rootcuatro = 16-19;       
constraint rootcinco = 21-48;       
constraint rootseis = 21-26;       

calibrate root=Offsetlognormal(21,25,4);
calibrate rootdos = Offsetlognormal(13,13.3,0.3); 
calibrate roottrees = Offsetlognormal(12,13,0.2);
calibrate rootcuatro = Offsetlognormal(12.8,13,0.6); 
calibrate rootcinco = Offsetlognormal(13.4,13.6,0.2);
calibrate rootseis = Offsetlognormal(10,11,0.8);
prset nodeagepr=calibrated;

prset topologypr=constraints(root, rootdos, roottrees, rootcuatro, rootcinco, rootseis);

mcmc ngen=30000000 nchains=4 printfreq=1000 samplefreq=1000 checkpoint=yes;

sump;
sumt;
```
