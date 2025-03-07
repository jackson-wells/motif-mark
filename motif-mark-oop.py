#!/usr/bin/python

import argparse
import os
import cairo
import re

COLORS = [
    (31/255, 119/255, 180/255, 0.75),  # Blue
    (255/255, 127/255, 14/255, 0.75),  # Orange
    (44/255, 160/255, 44/255, 0.75),   # Green
    (214/255, 39/255, 40/255, 0.75),   # Red
    (148/255, 103/255, 189/255, 0.75), # Purple
    (140/255, 86/255, 75/255, 0.75),   # Brown
    (227/255, 119/255, 194/255, 0.75), # Pink
    (188/255, 189/255, 34/255, 0.75),  # Yellow-green
    (23/255, 190/255, 207/255, 0.75)   # Cyan
]

class File:
	def __init__(self, file):
		'''Stores basic file info, to be used to for any file type taken in by the script'''
		self.file = file
		self.name = self._setName()
		self.extension = self._setExtension()
		self.path = self._setPath()
		
	def _setName(self):
		return os.path.splitext(self.file)[0]
	
	def _setExtension(self):
		return os.path.splitext(self.file)[1]
	
	def _setPath(self):
		return os.path.dirname(self.file)

class Fasta(File):
	def __init__(self, file):
		'''Inherits from file class, contains a list of objects that will contain all FASTA file information'''
		super().__init__(file)
		self.sequences = self._setSequences()

	def getSequenceCount(self):
		return len(self.sequences)

	def _setSequences(self):
		tempList = list()
		tempObject = None
		with open(self.file,"r") as fh:
			for line in fh:
				if not line.find(">"):  
					if tempObject is not None:
						tempList.append(tempObject)
					tempObject = Sequence("","")
					tempObject.setHeader(line.strip('\n'))
					tempObject.setGene()
				else:     
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
		tempList = list()
		with open(self.file,"r") as fh:
			for line in fh:
				tempList.append(line.strip('\n'))
		return tempList
	
	def getMotifs(self):
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
		temp = self.header.split(' ')[0]
		self.gene = temp[1:]

	def setBody(self, body):
		self.body = body
	
	def buildBody(self, sequence):
		self.body += sequence
	
	def setHeader(self, header):
		self.header = header

	def findExons(self):
		'''Identifies and stores any exons found within provided sequence'''
		tempList = list()
		temp = ""
		tempStart = 0
		for i, character in enumerate(self.body):
			if character.isupper():
				if temp == "":
					tempStart = i
				temp += character
			elif temp:
				tempList.append(Atrribute(temp,len(temp),tempStart,i-1))
				temp = ""
		if temp:
			tempList.append(Atrribute(temp,len(temp),tempStart,i-1))
		return tempList
	
	def findIntrons(self):
		'''Identifies and stores any introns found within provided sequence'''
		tempList = list()
		temp = ""
		tempStart = 0
		for i, character in enumerate(self.body):
			if character.islower():
				if temp == "":
					tempStart = i
				temp += character
			elif temp:
				tempList.append(Atrribute(temp,len(temp),tempStart,i-1))
				temp = ""

		if temp:
			tempList.append(Atrribute(temp,len(temp),tempStart,i-1))
		return tempList
	
	def findMotifs(self, motifs):
		'''Identifies and stores any motifs, from the input motif file, found within provided sequence'''
		tempSequence = self.body.lower().replace('u', 't')
		tempList = list()
		for motif in motifs:
			startPositions = list()
			stopPositions = list()
			tempCount = 0
			tempMotif = motif.lower().replace('u', 't')
			tempMotif = tempMotif.replace('y', '[tc]')
			pattern = re.compile(tempMotif)
			for match in re.finditer(pattern, tempSequence):
				startPositions.append(match.start())
				stopPositions.append(match.start() + len(motif) - 1)
				tempCount += 1
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
		self.occurences = temp
		return

class Draw:
	def __init__(self, width, height):
		self.width = width
		self.height = height
		self.surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width, height)
		self.ctx = cairo.Context(self.surface)

	def createCanvas(self):
		self.ctx.set_source_rgb(1, 1, 1)  
		self.ctx.paint()

	def drawImageTitle(self):
		# # Draw a title
		self.ctx.set_source_rgb(0, 0, 0)  # Black text
		self.ctx.select_font_face("Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
		self.ctx.set_font_size(32)

		# # Get text extenats for centering
		title = "Sequence Annotations"
		extents = self.ctx.text_extents(title)
		self.x = (self.width - extents.width) / 2  # Center horizontally
		self.ctx.move_to(self.x, 40)
		self.ctx.show_text(title) 

	def drawLegend(self, input_motifs):
		for i,motif in enumerate(input_motifs.motifs):
			self.y = 30 + i * 25  # Adjust vertical position for each entry

			# Draw color square
			self.ctx.set_source_rgba(*COLORS[i])
			self.ctx.rectangle(1000, self.y, 20, 20)
			self.ctx.fill()

			# Draw label text
			self.ctx.set_source_rgb(0, 0, 0)  # Black text
			self.ctx.set_font_size(16)
			self.ctx.move_to(1000 + 25, self.y + 20 - 5)  # Align text with square
			self.ctx.show_text(motif)

	def drawSequenceTitle(self, sequence):
		self.ctx.move_to(self.x, self.y)
		self.ctx.set_source_rgb(0, 0, 0)  # Black text
		self.ctx.select_font_face("Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
		self.ctx.set_font_size(24)
		tempTitle = "Gene: " + sequence.gene
		self.ctx.show_text(tempTitle)
		self.y += 50
		self.ctx.move_to(self.x, self.y)

	def align(self):
		self.x = 50
		self.y = 150

	def nextAnnotation(self):
		self.y += 100

	def drawIntrons(self, sequence):
		for i,intron in enumerate(sequence.introns):
			self.ctx.set_source_rgb(90/255, 90/255, 90/255) 
			self.ctx.line_to(self.x + intron.start + intron.length, self.y)
			self.ctx.stroke()
			self.ctx.move_to(self.x + intron.start + intron.length + sequence.exons[i-1].length, self.y)

	def drawExons(self, sequence):
		for exon in sequence.exons:
			# ctx.line_to(intron.start, intron.stop)
			self.ctx.rectangle(self.x + exon.start, self.y-10, exon.length, 20)
			self.ctx.set_source_rgb(90/255, 90/255, 90/255) 
			self.ctx.fill()
			self.ctx.stroke()
			self.ctx.move_to(50, 50)

	def drawMotifs(self, sequence):
		for i,motif in enumerate(sequence.motifs):
			for start in motif.start:
				self.ctx.rectangle(self.x + start, self.y-10, motif.length, 20)
				self.ctx.set_source_rgba(*COLORS[i]) 
				self.ctx.fill()
				self.ctx.stroke()
				self.ctx.move_to(50, 50)
	
	def outputImage(self, fileName):
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

	image = Draw(1200, (175*input_fasta.getSequenceCount()))

	image.createCanvas()
	image.drawImageTitle()
	image.drawLegend(input_motifs)
	image.align()

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

	#get motifs from file
	input_motifs = Motif(args.motif) # type: ignore

	#get FASTA sequences and headers from file
	input_fasta = Fasta(args.file) # type: ignore

	#Find attributes (exons, introns, motifs) in each FASTA sequence 
	findAttributes(input_fasta, input_motifs)

	#draw sequeneces with attributes 
	annotate(input_fasta, input_motifs)

if __name__ == "__main__":
    main()