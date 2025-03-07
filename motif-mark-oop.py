#!/usr/bin/python

import argparse
import os
import cairo
import re

#Constant list of colors 
COLORS = [
	#Blue
    (31/255, 119/255, 180/255, 0.75),  

	#Orange
    (255/255, 127/255, 14/255, 0.75),

	#Green  
    (44/255, 160/255, 44/255, 0.75),

	#Red   
    (214/255, 39/255, 40/255, 0.75), 

	#Purple  
    (148/255, 103/255, 189/255, 0.75), 

	#Brown
    (140/255, 86/255, 75/255, 0.75),   

	#Pink
    (227/255, 119/255, 194/255, 0.75), 

	#Puke green
    (188/255, 189/255, 34/255, 0.75),  

	#Cyan
    (23/255, 190/255, 207/255, 0.75)  
]

class File:
	def __init__(self, file):
		'''Stores basic file info, to be used to for any file type taken in by the script'''
		self.file = file
		self.name = self._setName()
		self.extension = self._setExtension()
		self.path = self._setPath()
		
	def _setName(self):
		'''Extract file name from provided input file'''
		return os.path.splitext(self.file)[0]
	
	def _setExtension(self):
		'''Extract file extension from provided input file'''
		return os.path.splitext(self.file)[1]
	
	def _setPath(self):
		'''Extract file path from provided input file'''
		return os.path.dirname(self.file)

class Fasta(File):
	def __init__(self, file):
		'''Inherits from file class, contains a list of objects that will contain all FASTA file information'''
		super().__init__(file)
		self.sequences = self._setSequences()

	def getSequenceCount(self):
		'''Gets number of sequences in the FASTA file'''
		return len(self.sequences)

	def _setSequences(self):
		'''Reads FASTA file and adds contents to list of sequence class objects'''
		tempList = list()
		tempObject = None
		with open(self.file,"r") as fh:
			for line in fh:
				#If line is the header line
				if not line.find(">"):  
					#If this a new header line, add object to list
					if tempObject is not None:
						tempList.append(tempObject)
					tempObject = Sequence("","")
					tempObject.setHeader(line.strip('\n'))
					tempObject.setGene()
				#If line is not header line
				else:
					#Check that object exists (catches invalid FASTA formatting (kind of))
					if tempObject is not None:                      
						tempObject.buildBody(line.strip('\n'))
			if tempObject is not None:
				tempList.append(tempObject)
		return tempList

class Motif(File):
	def __init__(self, file):
		'''Inherits from file class, contains a list of motifs found within the motif file'''
		super().__init__(file)
		self.motifs = self._setMotifs()

	def _setMotifs(self):
		'''Reads motif file and adds contents to list of strings'''
		tempList = list()
		with open(self.file,"r") as fh:
			for line in fh:
				tempList.append(line.strip('\n'))
		return tempList
	
	def getMotifs(self):
		'''returns list of motifs'''
		return self.motifs

class Sequence:
	def __init__(self, header, body):
		'''Contains all information within a FASTA file and empty lists to later contain sequence specific attributes'''
		self.header = header
		self.gene = ""
		self.body = body
		self.exons = list()
		self.introns = list()
		self.motifs = list()

	def setGene(self):
		'''Extracts gene name from sequence header'''
		temp = self.header.split(' ')[0]
		self.gene = temp[1:]

	def setBody(self, body):
		'''Updates sequence body according to single string'''
		self.body = body
	
	def buildBody(self, sequence):
		'''Updates sequence body by concatentation'''
		self.body += sequence
	
	def setHeader(self, header):
		'''Updates sequence header'''
		self.header = header

	def findExons(self):
		'''Identifies and stores any exons found within provided sequence'''
		tempList = list()
		temp = ""
		tempStart = 0
		#Loop over every character in sequence
		for i, character in enumerate(self.body):
			#If an Exon is detected
			if character.isupper():
				#Start of an exon
				if temp == "":
					#record start location
					tempStart = i
				#Build exon string
				temp += character
			#If lowercase/ end of exon
			elif temp:
				#add exon to list
				tempList.append(Atrribute(temp,len(temp),tempStart,tempStart+len(temp)))
				#reset exon string
				temp = ""
		if temp:
			tempList.append(Atrribute(temp,len(temp),tempStart,tempStart+len(temp)))
		return tempList
	
	def findIntrons(self):
		'''Identifies and stores any introns found within provided sequence'''
		tempList = list()
		temp = ""
		tempStart = 0
		#Loop over every character in sequence
		for i, character in enumerate(self.body):
			#if intron is detected
			if character.islower():
				#start of intron
				if temp == "":
					#record start location
					tempStart = i
				#build intron string
				temp += character
			#If uppercase/ end of intron
			elif temp:
				#add intron to list
				tempList.append(Atrribute(temp,len(temp),tempStart,tempStart+len(temp)))
				#reset intron string
				temp = ""
		if temp:
			tempList.append(Atrribute(temp,len(temp),tempStart,tempStart+len(temp)))
		return tempList
	
	def findMotifs(self, motifs):
		'''Identifies and stores any motifs, from the input motif file, found within provided sequence'''
		#create temporary lowercase sequence, replacing any U with T
		tempSequence = self.body.lower().replace('u', 't')
		tempList = list()
		#Loop over input motifs 
		for motif in motifs:
			startPositions = list()
			stopPositions = list()
			tempCount = 0
			#create temporary lowercase motif, replacing any U with T
			tempMotif = motif.lower().replace('u', 't')
			#Update any Y's to regex for T or C
			tempMotif = tempMotif.replace('y', '[tc]')
			pattern = re.compile(tempMotif)
			#Loop over identified occurences of motif in sequence
			for match in re.finditer(pattern, tempSequence):
				#Add start/stop positions to temp lists
				startPositions.append(match.start())
				stopPositions.append(match.start() + len(motif))
				tempCount += 1
			#If motif was found in sequence 
			if len(startPositions) > 0:
				tempAttribute = Atrribute(motif,len(motif),startPositions,stopPositions)
				tempAttribute.setOccurrences(tempCount)
				tempList.append(tempAttribute)
		return tempList
		
class Atrribute:
	def __init__(self, sequence, length, start, stop):
		'''Contains information for sequence attributes (exons, introns, motifs)'''
		self.sequence = sequence
		self.length = length
		self.start = start
		self.stop = stop
		self.occurences = 0

	def setOccurrences(self, temp):
		'''Updates occurence count (for motifs)'''
		self.occurences = temp
		
class Draw:
	'''Contains necessary components for drawing annotated sequences'''
	def __init__(self, width, height):
		self.width = width
		self.height = height
		self.surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width, height)
		self.ctx = cairo.Context(self.surface)

	def createCanvas(self):
		'''creates image canvas'''
		self.ctx.set_source_rgb(1, 1, 1)  
		self.ctx.paint()

	def drawImageTitle(self):
		'''Adds main title to image'''
		self.ctx.set_source_rgb(0, 0, 0)
		self.ctx.select_font_face("Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
		self.ctx.set_font_size(32)

		#Get text extenats for centering
		title = "Sequence Annotations"
		extents = self.ctx.text_extents(title)
		#Center title horizontally
		self.x = (self.width - extents.width) / 2  
		self.ctx.move_to(self.x, 40)
		self.ctx.show_text(title) 

	def drawLegend(self, input_motifs):
		'''Adds legend to image according to provided motifs'''
		for i,motif in enumerate(input_motifs.motifs):
			#move curson for each motif in legend
			self.y = 30 + i * 25  

			#Add color square corresponding to motif
			self.ctx.set_source_rgba(*COLORS[i])
			self.ctx.rectangle(1000, self.y, 20, 20)
			self.ctx.fill()

			#Add motif sequence 
			self.ctx.set_source_rgb(0, 0, 0) 
			self.ctx.set_font_size(16)
			#Align text with square
			self.ctx.move_to(1000 + 25, self.y + 20 - 5)  
			self.ctx.show_text(motif)

	def drawSequenceTitle(self, sequence):
		'''Adds sequnce title according to gene name in FASTA file'''
		self.ctx.move_to(self.x, self.y)
		self.ctx.set_source_rgb(0, 0, 0) 
		self.ctx.select_font_face("Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
		self.ctx.set_font_size(24)
		tempTitle = "Gene: " + sequence.gene
		self.ctx.show_text(tempTitle)
		#move down for sequence drawing 
		self.y += 50
		self.ctx.move_to(self.x, self.y)

	def alignMargin(self):
		'''aligns curson to write to draw sequences'''
		self.x = 50
		self.y = 150

	def nextAnnotation(self):
		'''moves curson down to draw next sequence'''
		self.y += 100

	def drawIntrons(self, sequence):
		'''Add introns by location on provide sequence'''
		for i,intron in enumerate(sequence.introns):
			self.ctx.set_source_rgb(90/255, 90/255, 90/255) #dark grey
			self.ctx.line_to(self.x + intron.start + intron.length, self.y)
			self.ctx.stroke()
			#move curson to after the exon that follows the current intron
			self.ctx.move_to(self.x + intron.start + intron.length + sequence.exons[i-1].length, self.y)

	def drawExons(self, sequence):
		'''Add exons by location on provide sequence'''
		#loop over identified exons in given sequence 
		for exon in sequence.exons:
			self.ctx.rectangle(self.x + exon.start, self.y-10, exon.length, 20)
			self.ctx.set_source_rgb(90/255, 90/255, 90/255) #dark grey
			self.ctx.fill()
			self.ctx.stroke()
			#reset curson to start of sequence
			self.ctx.move_to(50, 50)

	def drawMotifs(self, sequence):
		'''Annotate identifed motifs on provided sequence'''
		#Loop over motifs identifed to be present in sequence
		for i,motif in enumerate(sequence.motifs):
			#Loop over occurences of given motif
			for start in motif.start:
				#annotate motif using start position 
				self.ctx.rectangle(self.x + start, self.y-10, motif.length, 20)
				#Iteratively change color to delineate motifs
				self.ctx.set_source_rgba(*COLORS[i]) 
				self.ctx.fill()
				self.ctx.stroke()
				#reset curson to start of sequence
				self.ctx.move_to(50, 50)
	
	def outputImage(self, fileName):
		'''Output image according to provided filename'''
		self.surface.write_to_png(fileName)

def get_args():
    parser = argparse.ArgumentParser(description="python ./motif-mark-oop.py -f <FASTA file> -m <motif file> >")
    parser.add_argument("-f", "--file",help="Input FASTA file", type=str, required=True)
    parser.add_argument("-m", "--motif", help="valid motifs file", type=str, required=True)
    return parser.parse_args()

def findAttributes(input_fasta, input_motifs): 
	'''Using Fasta and Motif class objects, identifies introns, exons and motifs within each sequence form the input FASTA file'''
	for sequence in input_fasta.sequences:
		sequence.introns = sequence.findIntrons()
		sequence.exons = sequence.findExons()
		sequence.motifs = sequence.findMotifs(input_motifs.getMotifs())
	return

def annotate(input_fasta, input_motifs):
	'''Creates and outputs image containing annotated FASTA sequences'''
	#width allows space for legend to not overlap with sequences
	#height scales according to number of sequences in input FASTA file
	image = Draw(1200, (175*input_fasta.getSequenceCount()))
	image.createCanvas()

	image.drawImageTitle()
	image.drawLegend(input_motifs)
	image.alignMargin()

	#Loop over annotated sequences 
	for sequence in input_fasta.sequences:
		image.drawSequenceTitle(sequence)
		image.drawIntrons(sequence)
		image.drawExons(sequence)
		image.drawMotifs(sequence)
		image.nextAnnotation()

	image.outputImage(input_fasta.name + ".png")

	return

def main():
	args = get_args()

	#Get motifs from file
	input_motifs = Motif(args.motif)

	#Get FASTA sequences and headers from file
	input_fasta = Fasta(args.file) 

	#Find attributes (exons, introns, motifs) in each FASTA sequence 
	findAttributes(input_fasta, input_motifs)

	#Draw sequeneces with attributes 
	annotate(input_fasta, input_motifs)

if __name__ == "__main__":
    main()