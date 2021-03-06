<?xml version="1.0" encoding="utf-8"?>

<!-- $Id$ -->

<!--#########################################################################-->
<!--                                                                         -->
<!-- Project     Fittino - A SUSY Parameter Fitting Package                  -->
<!--                                                                         -->
<!-- File        Tool.xsl                                                    -->
<!--                                                                         -->
<!-- Description This file contains the xsl templates which specify how to   -->
<!--             display the configuration information of the used analysis  -->
<!--             tool                                                        -->
<!--                                                                         -->
<!-- Authors     Mathias Uhlenbrock  <uhlenbrock@physik.uni-bonn.de>         -->
<!--                                                                         -->
<!-- Licence     This program is free software; you can redistribute it      -->
<!--             and/or modify it under the terms of the GNU General Public  -->
<!--             License as published by the Free Software Foundation;       -->
<!--             either version 3 of he License, or (at your option) any     -->
<!--             later version.                                              -->
<!--                                                                         -->
<!--#########################################################################-->

<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:template match="Tool">
    <!-- Prints a headline with the tool's name and the individual configuration items -->
    <p>
      <!-- Headline -->
      <table>
        <tr><td><b>Tool</b></td><xsl:for-each select="*"><td><xsl:value-of select="local-name()"/></td></xsl:for-each></tr>
      </table>
    </p>
    <!-- Individual configuration items -->
    <xsl:apply-templates select="ContourPlotter |
                                 GeneticAlgorithmOptimizer |
                                 MarkovChainSampler |
                                 MinuitOptimizer |
                                 ParticleSwarmOptimizer |
                                 ScatterPlotter |
                                 SimpleSampler |
                                 SimulatedAnnealingOptimizer |
                                 SummaryPlotter |
                                 TreeSampler"/>
  </xsl:template>
  
  <!-- Individual configuration items of the tools -->
  
  <xsl:template match="ContourPlotter |
                       GeneticAlgorithmOptimizer |
                       MarkovChainSampler |
                       MinuitOptimizer |
                       ParticleSwarmOptimizer |
                       ScatterPlotter |
                       SimpleSampler |
                       SimulatedAnnealingOptimizer |
                       SummaryPlotter |
                       TreeSampler">
    <!-- Prints a list of the individual configuration items -->
    <p>
      <table>
        <xsl:for-each select="*">
          <!-- Create a new row for every configuration item -->
          <!-- TODO: Leading empty cells still have to be inserted manually -->
          <tr>
            <td></td><td><xsl:value-of select="local-name()"/></td><td class="cell-value"><xsl:value-of select="."/></td>
          </tr>
        </xsl:for-each>
      </table>
    </p>
  </xsl:template>

</xsl:stylesheet>
