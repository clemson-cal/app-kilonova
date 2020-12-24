use std::time::Instant;
use std::collections::HashMap;
use serde::{Serialize, Deserialize};




// ============================================================================
#[derive(Clone, Serialize, Deserialize)]
pub struct RecurringTask
{
    count: usize,
    next_time: f64,
}

#[derive(Clone, Serialize, Deserialize)]
pub struct Tasks
{
    pub write_checkpoint: RecurringTask,
    pub report_progress:  RecurringTask,
}

struct RunMonitor {    
    pub run_initiated:        Instant,
    pub last_report_progress: Instant,
    pub tasks_last_performed: Instant,
    pub call_count_this_run:  usize,
}




// ============================================================================
impl RecurringTask
{
    fn new() -> Self
    {
        Self{
            count: 0,
            next_time: 0.0,
        }
    }
    fn advance(&mut self, interval: f64)
    {
        self.count += 1;
        self.next_time += interval;
    }
}




// ============================================================================
impl Tasks
{
    fn new() -> Self
    {
        Self{
            write_checkpoint:     RecurringTask::new(),
            report_progress:      RecurringTask::new(),

        }
    }
}




// ============================================================================
impl RunMonitor {
    fn new() -> Self {
        Self{
            run_initiated:        Instant::now(),
            last_report_progress: Instant::now(),
            tasks_last_performed: Instant::now(),
            call_count_this_run:  0,
        }
    }
}
