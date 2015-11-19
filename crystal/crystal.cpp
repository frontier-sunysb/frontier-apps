/***************************************************************
FronTier is a set of libraries that implements different types of 
Front Traking algorithms. Front Tracking is a numerical method for 
the solution of partial differential equations whose solutions have 
discontinuities.  

Copyright (C) 1999 by The University at Stony Brook. 

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
****************************************************************/

/*
*				crystal.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*/

#include "crystal.h"

#define		MAX_NUM_VERTEX_IN_CELL		20  // ??

	/*  Local Application Function Declarations */

static void 	crystal_driver( Front*, C_CARTESIAN& );

char *in_name, *restart_state_name, *restart_name;
boolean RestartRun;
boolean ReadFromInput;
int RestartStep;

/*********************************************************************************
* Main function starts here
*********************************************************************************/
int main( int argc, char **argv )
{
	static Front front;            // main FronTier struct
	static F_BASIC_DATA f_basic;   // basic data in the main FronTier struct
	static LEVEL_FUNC_PACK level_func_pack; // level set initialization surface
	static VELO_FUNC_PACK velo_func_pack;  // velocity of the level set

	static CRT_PARAMS cRparams; // crystal parameters

	C_CARTESIAN       c_cartesian(front); // in crystal class

	FT_Init( argc, argv, &f_basic );  //  initialize f_basic 

	/* Initialize basic computational data */

	f_basic.size_of_intfc_state = sizeof(STATE); // init state of app
	
	/*Initialize Petsc */
	PetscInitialize( &argc, &argv, PETSC_NULL, PETSC_NULL);

        // get some data from f_basic
	in_name      		= f_basic.in_name;   // input filename
	restart_state_name      = f_basic.restart_state_name; // 
        restart_name 		= f_basic.restart_name; // restart directory
        RestartRun   		= f_basic.RestartRun;  // restart or not boolean
        ReadFromInput   	= f_basic.ReadFromInput;
	RestartStep 		= f_basic.RestartStep;

	sprintf( restart_state_name, "%s/state.ts%s", restart_name,
			right_flush(RestartStep,7)); // create string restart_state_name

	sprintf( restart_name, "%s/intfc-ts%s", restart_name,
			right_flush(RestartStep,7)); // create string restart_name

#if defined(__MPI__)
	if (pp_numnodes() > 1) // FronTier wrapper for MPI
	{
            sprintf( restart_name,"%s-nd%s",restart_name,
				right_flush(pp_mynode(),4));
            sprintf( restart_state_name,"%s-nd%s",restart_state_name,
				right_flush(pp_mynode(),4));
	}
#endif /* defined(__MPI__) */

	if (!ReadFromInput)  // need input file
	{
	    (void) printf("ERROR: Input file needed!\n");
	    clean_up(ERROR);
	}

	FT_ReadSpaceDomain( in_name, &f_basic ); // setup domain

	FT_StartUp( &front, &f_basic ); // 

	FT_InitDebug( in_name ); // define debug variables

	cRparams.dim = f_basic.dim;        // crystal
	front.extra2 = (POINTER)&cRparams; // FronTier accesses crystal params 

	if (!RestartRun)  // if first run
	{
	    /* Initialize interface through level function */
	    setInitialIntfc( &front, &level_func_pack, in_name );
	    FT_InitIntfc( &front, &level_func_pack);
	}
	else
	    read_restart_params( f_basic.dim, in_name, &front);


	read_crt_dirichlet_bdry_data( in_name, &front, f_basic);


	if (f_basic.dim == 2)

	    FT_ClipIntfcToSubdomain( &front ); // MPI domain decomposition

	if ( debugging( "init" ) ) // FronTier function; add "init" to debug options in the input file for running the application
	{
	    if (f_basic.dim == 2)
	    {
		char xg_name[100];  // for xgraph plotting 2D
		sprintf( xg_name, "init_intfc-%d", pp_mynode() );
		xgraph_2d_intfc( xg_name, front.interf );
	    }
	    else if ( f_basic.dim == 3 )
	    {
		char dname[100];  // for geomview 3D viz
		sprintf( dname, "init_intfc-%s", right_flush(pp_mynode(),4) );
		gview_plot_color_interface(dname,front.interf,YES);
	    }
	}

	read_crystal_params( in_name, &cRparams ); // x-tal non-member function init params

	FT_ReadTimeControl( in_name, &front ); // read time interval for plotting/print

	/* Initialize velocity field function */

	velo_func_pack.func_params = (POINTER)&cRparams; // FronTier uses this
	velo_func_pack.func = NULL;

	FT_InitFrontVeloFunc( &front, &velo_func_pack ); // FronTier init of intfc velo

        c_cartesian.initMesh(); // crystal init

        c_cartesian.initMovieVariables(); // what variable to do viz on

	if ( debugging("sample_solute") )
            c_cartesian.initSampleSolute( in_name );

        /*---------------------------------------
         key state update definition for FronTier
         ---------------------------------------*/
        front._point_propagate = crystal_point_propagate; 

 
	if (RestartRun)  // read state from input files
	{
	    c_cartesian.readFrontInteriorStates( restart_state_name );
	}
	else // init the states 
	{
	    initFrontStates( &front );
	    c_cartesian.setInitialCondition();
	}
	cRparams.field->vel = NULL; 	// No convection

        /*---------------------------------------
	 Propagate the front 
         ---------------------------------------*/
	crystal_driver( &front, c_cartesian );

	PetscFinalize();

	clean_up(0);
}
/*********************************************************************************
* Main function ends here
*********************************************************************************/

static  void crystal_driver(
        Front *front,
	C_CARTESIAN &c_cartesian )
{
        double CFL;

        // choose FronTier redistribution function
	Curve_redistribution_function(front) = expansion_redistribute;

        CFL = Time_step_factor(front);

        // draw initial state
	c_cartesian.crystalDraw();

        if (!RestartRun)
        {
	    FT_ResetTime(front);
	    FT_Propagate(front);

	    c_cartesian.solve( front->dt );

	    c_cartesian.timeStepAnalysis( NO );

	    FT_SetTimeStep( front );

	    c_cartesian.setAdvectionDt();

	    front->dt = std::min( front->dt, CFL*c_cartesian.max_dt );

            // save states to file 
            FT_Save( front );

	    c_cartesian.printFrontInteriorStates();

	    FT_SetOutputCounter( front );

            // save data for movie
            FT_Draw( front );
        }
        else // if it is a restart
        {
	    FT_SetOutputCounter( front );
        }

        // choose dt among dt's for print, propagate, movie, solver
	FT_TimeControlFilter( front );

	FT_PrintTimeStamp( front );

        for (;;)
        {
	    FT_Propagate(front);
	    c_cartesian.solve(front->dt);
	    FT_AddTimeStepToCounter(front);

	    FT_SetTimeStep(front);
	    front->dt = std::min(front->dt,CFL*c_cartesian.max_dt);

            if (FT_IsSaveTime(front))
	    {
                FT_Save(front);
		c_cartesian.printFrontInteriorStates();
	    }
            if (FT_IsDrawTime(front))
	    {
                FT_Draw(front);
                c_cartesian.crystalDraw();
	    }
            if (FT_TimeLimitReached(front) || bdryReached(front))
	    {
	    	c_cartesian.timeStepAnalysis(YES);
	    	FT_PrintTimeStamp(front);
                FT_Draw(front);
                FT_Save(front);
		c_cartesian.printFrontInteriorStates();
                break;
	    }
            if (FT_IsDrawTime(front))
	    	c_cartesian.timeStepAnalysis(YES);
	    else
	    	c_cartesian.timeStepAnalysis(NO);
	    FT_TimeControlFilter(front);
	    FT_PrintTimeStamp(front);
        }
}       /* end crystal_driver */
