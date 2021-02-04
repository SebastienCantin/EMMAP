// STAR-CCM+ macro: Init_Coupling.java.java
// Written by STAR-CCM+ 12.04.010
// At the second temporal iteration, the coupling is added as a sink term in the Region considered.
package macro;

import java.util.*;

import star.common.*;
import star.base.neo.*;
import star.species.*;

public class Init_Coupling.java extends StarMacro {

  public void execute() {
    execute0();
  }

  private void execute0() {

    Simulation simulation_0 = 
      getActiveSimulation();

    simulation_0.getSimulationIterator().step(1);

     Region region_0 = 
      simulation_0.getRegionManager().getRegion("FLUID"); /* Default name of Region is "FLUID" */

    SpeciesUserSource speciesUserSource_0 = 
      region_0.getValues().get(SpeciesUserSource.class);

    ScalarProfile scalarProfile_0 = 
      speciesUserSource_0.getMethod(CompositeArrayProfileMethod.class).getProfile(1);

    XyzInternalTable xyzInternalTable_0 = 
      ((XyzInternalTable) simulation_0.getTableManager().getTable("RateCondensation")); /* Create a table RateCondensation in STAR-CCM+ */

    scalarProfile_0.getMethod(XyzTabularScalarProfileMethod.class).setTable(xyzInternalTable_0);

    scalarProfile_0.getMethod(XyzTabularScalarProfileMethod.class).setData("User SumCondensation");
 
    simulation_0.getSimulationIterator().run();
    
   
    
  }
}
