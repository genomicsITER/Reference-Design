package org.broadinstitute.gatk.queue.qscripts.examples

import org.broadinstitute.gatk.queue.QScript
import org.broadinstitute.gatk.queue.extensions.gatk._
import org.broadinstitute.gatk.utils.commandline.Hidden
import org.broadinstitute.gatk.utils.commandline._
import org.broadinstitute.gatk.queue.util.QScriptUtils
import org.broadinstitute.gatk.queue.function.ListWriterFunction
//import org.broadinstitute.gatk.utils.baq.BAQ.CalculationMode

/*
How to run:  java -jar Queue.jar -R Homo_sapiens_assembly19.fasta -I test_dupRem.bam -indels Homo_sapiens_assembly19.known_indels.vcf -S ExampleIndelRealigner.scala -l DEBUG -run -jobRunner CMPShell
*/
class ExampleIndelRealigner extends QScript {
  // Create an alias 'qscript' to be able to access variables
  // in the ExampleUnifiedGenotyper.
  // 'qscript' is now the same as 'ExampleUnifiedGenotyper.this'
  qscript =>


  // Required arguments.  All initialized to empty values.

  @Input(doc="The reference file for the bam files.", shortName="R")
  var referenceFile: File = _ // _ is scala shorthand for null

  @Input(doc="Bam file to indel realigner.", shortName="I")
  var bamFile: File = _


  //Optional
  @Input(doc="extra VCF files to use as reference indels for Indel Realignment", fullName="extra_indels", shortName="indels", required=false)
  var indels: Seq[File] = Seq()


  /****************************************************************************
  * Hidden Parameters
  ****************************************************************************/
  @Hidden
  @Argument(doc="How many ways to scatter/gather", fullName="scatter_gather", shortName="sg", required=false)
  var nContigs: Int = -1

  @Argument(doc="Define the default platform for Count Covariates -- useful for techdev purposes only.", fullName="default_platform", shortName="dp", required=false)
  var defaultPlatform: String = ""


  // This trait allows us set the variables below in one place,
  // and then reuse this trait on each CommandLineGATK function below.

  // General arguments to GATK walkers
  trait CommandLineGATKArgs extends CommandLineGATK {
    this.reference_sequence = qscript.referenceFile
  }

  case class clean (inBam: File, tIntervals: File, outBam: File) extends IndelRealigner with CommandLineGATKArgs {
    this.input_file :+= inBam
    if (qscript.indels != null)
       this.known ++= qscript.indels
    this.targetIntervals = tIntervals
    this.out = outBam
    this.filter_bases_not_stored = true
    this.scatterCount = qscript.nContigs
    this.memoryLimit = 4
  }

  
  def script() {
      nContigs = 16//QScriptUtils.getNumberOfContigs(bamFile)
      val dupRemBam = bamFile
      val targetIntervals  = swapExt(dupRemBam.getParent, dupRemBam, "_dupRem.bam", "_realn.intervals")
      val realignBam = swapExt(dupRemBam.getParent, dupRemBam, "_dupRem.bam", "_realign.bam")
      add(clean(dupRemBam, targetIntervals, realignBam))
    }
}
