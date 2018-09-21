The following codes are used for lumping bone biology model.
Basically Slin.m and Slum.m are the only codes to be edited

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%  Lumped version of
%    Peterson MC, Riggs MM (2010) Bone 46:49-63
%                        +
%    Peterson MC, Riggs MM (2012) CPT Pharmacometrics Syst Pharmacol 1:e14
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

Run file = RUN_lump.m

1. linearize the original nonlinear model

        Slin

	Model_parameter_values0
	IC_setting

	<Only for original nonlinear system>
	crit_unlumped_PD
	pkpdfun
	def_ode
	
	<Only for linearizing denosumab PK>
	Kmat_induc_pk

	<Only for linearizing bone system>
	Kmat_induc_PD

	ME_solution

2. lump the linearized model

        Slum

	Model_parameter_values_SE
	crit_unlumped_PD_SE
	pkpdfun_SE

	vector2matrix
	matrix2vector
	OBJV_function_K_BMD
	nearestNeighbour
	OBJV_function_vpc

        Criterion
