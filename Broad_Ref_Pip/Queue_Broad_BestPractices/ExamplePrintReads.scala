package org.broadinstitute.gatk.queue.qscripts.examples

import org.broadinstitute.gatk.queue.QScript
import org.broadinstitute.gatk.queue.extensions.gatk._
import org.broadinstitute.gatk.utils.commandline.Hidden
import org.broadinstitute.gatk.utils.commandline._
import org.broadinstitute.gatk.queue.util.QScriptUtils
import org.broadinstitute.gatk.queue.function.ListWriterFunction
import org.broadinstitute.gatk.utils.baq.BAQ.CalculationMode
/* How to run:
java -jar Queue.jar -R Homo_sapiens_assembly19.fasta -I test_realign.bam -B test_recal.grp -l DEBUG -S ExamplePrintReads.scala -run -jobRunner CMPShell
*/

class ExamplePrintReads extends QScript {
  // Create an alias 'qscript' to be able to access variables
  // in the ExampleUnifiedGenotyper.
  // 'qscript' is now the same as 'ExampleUnifiedGenotyper.this'
  qscript =>


  // Required arguments.  All initialized to empty values.

  @Input(doc="The reference file for the bam files.", shortName="R")
  var referenceFile: File = _ // _ is scala shorthand for null

  @Input(doc="Bam file to genotype.", shortName="I")
  var bamFile: File = _

  @Input(doc="BQSR recal file", fullName="BQSR", shortName="B", required=true)
  var recalFile: File = _

  @Input(doc="dbsnp ROD to use (must be in VCF format)", fullName="dbsnp", shortName="D", required=false)
   var dbSNP: Seq[File] = Seq()
    
  
  // The following arguments are all optional.


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


  
  case class recal (inBam: File, inRecalFile: File, outBam: File) extends PrintReads with CommandLineGATKArgs {
    this.input_file :+= inBam
    this.BQSR = inRecalFile
    this.out = outBam
    this.scatterCount = qscript.nContigs
    this.memoryLimit = 4
  }


  def script() {
      nContigs = 36//QScriptUtils.getNumberOfContigs(bamFile)
      val dedupedBam = qscript.bamFile
      val preRecalFile  = qscript.recalFile
      val recalBam = swapExt(dedupedBam.getParent, dedupedBam, "_realign.bam", "_final.bam")
      add(recal(dedupedBam, preRecalFile, recalBam))
    }
}
