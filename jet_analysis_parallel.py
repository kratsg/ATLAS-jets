from atlas_jets import *
filename = '/Users/kratsg/Dropbox/UChicagoSchool/DavidMiller/Data/PileupSkim_TTbar_14TeV_MU80.root'
directory = 'TTbar_14TeV_MU80/'
tree = 'mytree'

from multiprocessing import Process, Queue, cpu_count

class Consumer(Process):
  def __init__(self, tasks, results):
    Process.__init__(self)
    self.tasks = tasks
    self.results = results

  def run(self):
    proc_name = self.name
    while True:
      next_task = self.tasks.get()
      if next_task is None:
        # Poison pill means we should exit
        break
      result = next_task()
      self.results.put(result)
    return

class Task(object):
  def __init__(self, oEvent):
    self.oEvent = oEvent
  def __call__(self):
    grid = OfflineJets.Grid(cell_resolution=0.2)
    grid.add_event(self.oEvent)
    return self.oEvent
  def __str__(self):
    pass

# Must use this to avoid spawning multiple processes and screwing up name-scope etc...
#     see python docs style guide for particular reasons, I'm just bullshitting
#       but we really do need this
if __name__ == '__main__':

  oEvents = OfflineJets.Events(rootfile)
  processed_oEvents = []

  # Establish communication queues
  tasks, results = Queue(), Queue()

  # Build up the processes
  num_consumers = cpu_count()
  num_tasks = 0
  consumers = [ Consumer(tasks, results) for i in range(num_consumers) ]
  for oEvent in oEvents:
    tasks.put(Task(oEvent))
    num_tasks += 1

  # poison pill to stop processing
  for i in range(num_consumers):
    tasks.put(None)

  # start processing
  for c in consumers:
    c.start()

  # wait for all tasks to complete
  while num_tasks:
    processed_oEvents.append(results.get())#this should be blocking
    num_tasks -= 1

  # at this point, processed_events is the events array that contains calculated events with no issue
  analysis = Analysis(offline_events = processed_oEvents)
  analysis.Efficiency()
