package org.broadinstitute.gatk.queue.qscripts.examples

import org.broadinstitute.gatk.queue.QScript
import org.broadinstitute.gatk.queue.extensions.gatk._
import org.broadinstitute.gatk.utils.commandline.Hidden
import org.broadinstitute.gatk.utils.commandline._
import org.broadinstitute.gatk.queue.util.QScriptUtils
import org.broadinstitute.gatk.queue.function.ListWriterFunction
import org.broadinstitute.gatk.utils.variant.GATKVCFIndexType
import org.broadinstitute.gatk.tools.walkers.haplotypecaller.ReferenceConfidenceMode

//import org.broadinstitute.gatk.utils.baq.BAQ.CalculationMode

/*
How to run:  java -jar Queue.jar -R Homo_sapiens_assembly19.fasta -I test_final.bam -S ExampleHaplotypeCaller.scala -l DEBUG -run -jobRunner CMPShell
*/
class ExampleHaplotypeCaller extends QScript {
  // Create an alias 'qscript' to be able to access variables
  // in the ExampleUnifiedGenotyper.
  // 'qscript' is now the same as 'ExampleUnifiedGenotyper.this'
  qscript =>


  // Required arguments.  All initialized to empty values.

  @Input(doc="The reference file for the bam files.", shortName="R")
  var referenceFile: File = _ // _ is scala shorthand for null

  @Input(doc="Bam file to indel realigner.", shortName="I")
  var bamFile: File = _

  // The following arguments are all optional.

  @Input(doc="an intervals file to be used by GATK - output bams at intervals only", fullName="gatk_interval_file", shortName="intervals", required=false)
  var intervals: File = _


  /****************************************************************************
  * Hidden Parameters
  ****************************************************************************/
  @Hidden
  @Argument(doc="How many ways to scatter/gather", fullName="scatter_gather", shortName="sg", required=false)
  var nContigs: Int = -1

  // This trait allows us set the variables below in one place,
  // and then reuse this trait on each CommandLineGATK function below.

  // General arguments to GATK walkers
  trait CommandLineGATKArgs extends CommandLineGATK {
    this.reference_sequence = qscript.referenceFile
  }

  case class varcall (inBam: File, outVCF: File) extends HaplotypeCaller with CommandLineGATKArgs {
    this.input_file :+= inBam
    this.out = outVCF
    this.ip = 50
    this.emitRefConfidence = ReferenceConfidenceMode.GVCF
    this.variant_index_type =  GATKVCFIndexType.LINEAR
    this.variant_index_parameter = 128000
    this.intervals = if (qscript.intervals == null) Nil else List(qscript.intervals)
    this.scatterCount = qscript.nContigs
    this.memoryLimit = 4
  }

  
  def script() {
      nContigs = 16//QScriptUtils.getNumberOfContigs(bamFile)
      val recalBam = qscript.bamFile
      val finalVCF = swapExt(recalBam.getParent, recalBam, "_final.bam", "_HaplotypeCaller.g.vcf")
      add(varcall(recalBam, finalVCF))
    }
}
