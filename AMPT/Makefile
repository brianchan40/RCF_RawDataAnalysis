CENTRALITIES = 0 1 2 3 4 5 6 7 8
#CENTRALITIES = 8
#PERIOD_5 = 161
#PERIOD_4 = 1 2 3 4 5 6 7 8
PERIOD_4 = 1
#LAMBDA = lamrot
LAMBDA = antilam
#LUM = high mid low nor
LUM = low

submit_gamma:
	rm -rf ./pla_*.zip  ./pla_*.package
	rm -rf ./output/*
	for lum in $(LUM) ; do for lam in $(LAMBDA) ; do for cen in $(CENTRALITIES) ; do ./submit_gamma.sh $$cen 0 $$lam $$lum; done ; done ; done;

submit_v2:
	rm -rf ./pla_*.zip  ./pla_*.package
	rm -rf ./output/*
	for lum in $(LUM) ; do for lam in $(LAMBDA) ; do for cen in $(CENTRALITIES) ; do ./submit_v2.sh $$cen 0 $$lam $$lum; done ; done ; done;

analysis:
	for daynum in $(PERIOD_4); do \
		./run_gamma.csh 1 1 $$daynum 1234; \
		mv 1234cen1.weight_112_module_new.root cen1.weight_112_module.root; \
		mv cen1.v2_fullEP_eff_pT02_module.root cen1.v2_fullEP_eff_pT02_module_old.root; \
		./run_gamma.csh 1 1 $$daynum 1234; \
		mv cen1.weight_112_module.root cen1.weight_112_module_old.root; \
		mv 1234cen1.weight_112_module_new.root cen1.weight_112_module.root; \
		./run_gamma.csh 1 0 $$daynum 1234; \
		mv cen* ./Results/cen1/Data160_2/; \
	done;

prep:
	for cen in $(CENTRALITIES) ; do mkdir ./condor_files/Gamma_$${cen} ; done ;

combine:
	for lam in $(LAMBDA) ; do for cen in $(CENTRALITIES) ; do \
                hadd ./Results_$${lam}_18/cen$$cen.weight_112_module.root ./temp/*/*cen$$cen.weight_112_module_new.root; \
                hadd ./Results_$${lam}_18/cen$$cen.v2_fullEP_eff_pT02_module.root ./temp/*/*cen$$cen.v2_fullEP_eff_pT02_module.root; \
        done; done;

combine1:
	for lam in $(LAMBDA) ; do for cen in $(CENTRALITIES) ; do \
		mv ./Results_$${lam}_18/cen$$cen.weight_112_module.root ./Results_$${lam}_18/cen$$cen.weight_112_module_old.root; \
		mv ./Results_$${lam}_18/cen$$cen.v2_fullEP_eff_pT02_module.root ./Results_$${lam}_18/cen$$cen.v2_fullEP_eff_pT02_module_old.root; \
		hadd ./Results_$${lam}_18/cen$$cen.weight_112_module.root ./temp/*/*cen$$cen.weight_112_module_new.root; \
		hadd ./Results_$${lam}_18/cen$$cen.v2_fullEP_eff_pT02_module.root ./temp/*/*cen$$cen.v2_fullEP_eff_pT02_module.root; \
	done; done ;

reversecombine1:
	for lam in $(LAMBDA) ; do for cen in $(CENTRALITIES) ; do \
                mv ./Results_$${lam}_18/cen$$cen.weight_112_module_old.root ./Results_$${lam}_18/cen$$cen.weight_112_module.root; \
                mv ./Results_$${lam}_18/cen$$cen.v2_fullEP_eff_pT02_module_old.root ./Results_$${lam}_18/cen$$cen.v2_fullEP_eff_pT02_module.root; \
        done; done ;

combine2:
	for lam in $(LAMBDA) ; do for cen in $(CENTRALITIES) ; do \
                hadd ./Results_$${lam}_18/cen$$cen.gamma112_fullEP_eff_pT02_module.root ./temp/*/*cen$$cen.gamma112_fullEP_eff_pT02_module.root; \
        done; done ;

c:
	rm *.session.xml
	rm -rf condor_files/*
	rm sched*
	rm -r output/*/*_debug
	mkdir ./condor_files/Gamma112_0
	mkdir ./condor_files/Gamma112_1
	mkdir ./condor_files/Gamma112_2
	mkdir ./condor_files/Gamma112_3
	mkdir ./condor_files/Gamma112_4
	mkdir ./condor_files/Gamma112_5
	mkdir ./condor_files/Gamma112_6
	mkdir ./condor_files/Gamma112_7
	mkdir ./condor_files/Gamma112_8

dataquery:
	get_file_list.pl -keys 'path,filename' -cond 'production=P17ih,trgsetupname=AuAu54_production_2017,filetype=daq_reco_MuDst,filename~st_physics_181,storage=HPSS' -limit 500
kill:
	condor_rm brian40	
check:
	@condor_q brian40 | tail -10

.PHONY: clean submit_phi submit_ep dataquery kill check test test_loop

