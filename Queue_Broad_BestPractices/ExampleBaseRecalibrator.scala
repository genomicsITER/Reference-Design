package org.broadinstitute.gatk.queue.qscripts.examples

import org.broadinstitute.gatk.queue.QScript
import org.broadinstitute.gatk.queue.extensions.gatk._
import org.broadinstitute.gatk.utils.commandline.Hidden
import org.broadinstitute.gatk.utils.commandline._
import org.broadinstitute.gatk.queue.util.QScriptUtils
import org.broadinstitute.gatk.queue.function.ListWriterFunction

/*
How to run:  java -jar Queue.jar -R Homo_sapiens_assembly19.fasta -I test_realign.bam -D Homo_sapiens_assembly19.dbsnp.vcf -S ExampleBaseRecalibrator.scala -l DEBUG -run -jobRunner CMPShell
*/
class ExampleBaseRecalibrator extends QScript {
  // Create an alias 'qscript' to be able to access variables
  // in the ExampleUnifiedGenotyper.
  // 'qscript' is now the same as 'ExampleUnifiedGenotyper.this'
  qscript =>


  // Required arguments.  All initialized to empty values.

  @Input(doc="The reference file for the bam files.", shortName="R")
  var referenceFile: File = _ // _ is scala shorthand for null

  @Input(doc="Bam file to genotype.", shortName="I")
  var bamFile: File = _

  @Input(doc="dbsnp ROD to use (must be in VCF format)", fullName="dbsnp", shortName="D", required=true)
  var dbSNP: Seq[File] = Seq()

  
  // The following arguments are all optional.

  @Input(doc="an intervals file to be used by GATK - output bams at intervals only", fullName="gatk_interval_file", shortName="intervals", required=false)
  var intervals: File = _


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


  case class cov (inBam: File, outRecalFile: File) extends BaseRecalibrator with CommandLineGATKArgs {  
    this.input_file :+= inBam
    this.knownSites ++= qscript.dbSNP
    this.ip = 50
    this.out = outRecalFile
    if (!defaultPlatform.isEmpty) this.default_platform = defaultPlatform

    this.intervals = if (qscript.intervals == null) Nil else List(qscript.intervals)
    this.scatterCount = qscript.nContigs
      // Set the memory limit to 2 gigabytes on each command.
    this.memoryLimit = 4
  }

  def script() {
      nContigs = 16//QScriptUtils.getNumberOfContigs(bamFile)
      val dedupedBam = bamFile
      val preRecalFile  = swapExt(dedupedBam.getParent, dedupedBam, "_realign.bam", "_recal.grp")
      add(cov(dedupedBam, preRecalFile))
    }
}
