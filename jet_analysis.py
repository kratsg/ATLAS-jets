import root_jets as Jets
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
  def __init__(self, event):
    self.event = event
  def __call__(self):
    grid = Jets.Grid(pixel_resolution=0.02)
    grid.add_event(self.event)
    return self.event
  def __str__(self):
    pass

class JetAnalysis:
  '''A good holder for all the jet analysis we might do - computationally'''
  def __init__(self):
    pass

  def efficiency(self, events=[]):
    pass

# Must use this to avoid spawning multiple processes and screwing up name-scope etc...
#     see python docs style guide for particular reasons, I'm just bullshitting
#       but we really do need this
if __name__ == '__main__':
  # load up the events
  filename = '/Users/kratsg/Dropbox/UChicagoSchool/DavidMiller/Data/gFEXSlim_ttbar_zprime3000.root'
  events = Jets.Events(filename=filename)
  processed_events = []
  events.load()

  # Establish communication queues
  tasks, results = Queue(), Queue()

  # Build up the processes
  num_consumers = cpu_count()
  num_tasks = 0
  consumers = [ Consumer(tasks, results) for i in range(num_consumers) ]
  for event in events:
    tasks.put(Task(event))
    num_tasks += 1

  for jet in events[0]:
    print jet

  # poison pill to stop processing
  for i in range(num_consumers):
    tasks.put(None)

  # start processing
  for c in consumers:
    c.start()

  # wait for all tasks to complete
  while num_tasks:
    processed_events.append(results.get())#this should be blocking
    num_tasks -= 1

  # at this point, processed_events is the events array that contains calculated events with no issue
