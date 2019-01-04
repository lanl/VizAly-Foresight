#!/usr/bin/python

import sys

if py3:
    import io
else:
    import StringIO as io

dbConnLoaded = ""

try:
	import apsw
	dbConnLoaded = "apsw"
	print("using aspw module for sqlite")
except ImportError:
	try:
		import sqlite3
		dbConnLoaded = "sqlite3"
		print("using sqlite3 module for sqlite")
	except ImportError:
		print ("You need either apsw or sqlite3")
		sys.exit(0)


class GioSqlite3:

	def __init__(self):
		if dbConnLoaded == "sqlite3":
			self.conn = sqlite3.connect(":memory:")
		else:
			self.conn = apsw.Connection(":memory:")

			output=io.StringIO()

			self.shell = apsw.Shell(stdout=output, db=self.conn)
			self.shell.process_command(".mode csv")
			self.shell.process_command(".headers on")

	def __del__(self):
		self.closeConn()


	def loadGIOSqlite(self, sharedObjectPath):
		try:
			if dbConnLoaded == "sqlite3":
				self.conn.enable_load_extension(True)
				self.conn.load_extension(sharedObjectPath)
			else:
				self.conn.enableloadextension(True)
				self.conn.loadextension(sharedObjectPath)
		except Exception:
			print ("Could not load shared object", sharedObjectPath, "!")
			return -1


	def createTable(self, tableName, inputfile):
		query = "CREATE VIRTUAL TABLE " + tableName + " USING GenericIO('"+ inputfile +"')"  
		try:
			self.conn.cursor().execute(query)
		except Exception:
			print ("Could not create table", tableName, "!")
			return -1


	def runQueryInteractive(self, queryString):
		for row in self.conn.cursor().execute(queryString):
			print (row)


	def runQueryOutputFile(self, queryString, outputFilename):
		target = open(outputFilename, 'w')

		for row in self.conn.cursor().execute(queryString):
			target.write( str(row) + '\n')


	def runQueryOutputString(self, queryString): 		
		outputString = "" 				
		cursor = self.conn.cursor().execute(queryString) 		
		results = cursor.fetchall() 		
		for row in results: 			
			outputString = outputString + str( row ) + " "'\n'

		return outputString


	def runQueryOutputList(self, queryString):
		outputString = ""
		cursor = self.conn.cursor().execute(queryString)
		row = cursor.fetchone()

		resultsList = []
		while row is not None:
			resultsList.append(row[0])
			row = cursor.fetchone()

		return resultsList



	def getNumRanks(self, tableName):
		query = "SELECT MAX(_rank) FROM " + tableName
		cursor = self.conn.cursor().execute(query)
		return cursor.fetchone()[0]


	def closeConn(self):
		self.conn.close()