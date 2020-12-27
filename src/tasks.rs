use std::time::Instant;
use serde::{Serialize, Deserialize};




/**
 * A task, or side-effect, such as reporting, analysis, or data output
 */
#[derive(Clone, Serialize, Deserialize)]
pub struct RecurringTask {

    /// The number of times this task has been performed
    count: usize,

    /// The next simulation time at which this task is set to be performed
    next_time: f64,

    /// The last clock time when this task was performed
    #[serde(skip, default = "Instant::now")]
    last_performed: Instant,
}




/**
 * All the tasks that are used in this application
 */
#[derive(Clone, Serialize, Deserialize)]
pub struct Tasks {

    /// Write a snapshot of the full simulation
    pub write_checkpoint: RecurringTask,

    /// Output the primitive and geometric quantities for plotting and
    /// post-processing
    pub write_primitives: RecurringTask,

    /// Summarize the simulation performance
    pub report_progress:  RecurringTask,
}




// ============================================================================
impl RecurringTask
{

    /**
     * Create a fresh recurring task which is first due at t = 0.0.
     */
    pub fn new() -> Self {
        Self{
            count: 0,
            next_time: 0.0,
            last_performed: Instant::now(),
        }
    }

    /**
     * Mark the task as having just been performed, and schedule it to happen
     * again after the given time interval.
     */
    pub fn advance(&mut self, interval: f64) {
        self.count += 1;
        self.next_time += interval;
    }
}




// ============================================================================
impl Tasks
{
    pub fn new() -> Self {
        Self{
            write_checkpoint: RecurringTask::new(),
            write_primitives: RecurringTask::new(),
            report_progress: RecurringTask::new(),
        }
    }
}
