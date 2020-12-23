use std::time::Instant;
use std::collections::HashMap;
use clap::App;
use kind_config::Form;




// ============================================================================
#[derive(Clone, hdf5::H5Type)]
#[repr(C)]
pub struct RecurringTask
{
    count: usize,
    next_time: f64,
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
#[derive(Clone)]
pub struct Tasks
{
    pub write_checkpoint:     RecurringTask,
    pub report_progress:      RecurringTask,
    pub run_initiated:        Instant,
    pub last_report_progress: Instant,
    pub tasks_last_performed: Instant,
    pub call_count_this_run:  usize,
}




// ============================================================================
impl From<Tasks> for Vec<(String, RecurringTask)>
{
    fn from(tasks: Tasks) -> Self {
        vec![
            ("write_checkpoint".into(), tasks.write_checkpoint),
            ("report_progress".into(), tasks.report_progress),
        ]
    }
}

impl From<Vec<(String, RecurringTask)>> for Tasks
{
    fn from(a: Vec<(String, RecurringTask)>) -> Tasks {
        let task_map: HashMap<_, _> = a.into_iter().collect();
        let mut tasks = Tasks::new();
        tasks.write_checkpoint   = task_map.get("write_checkpoint")  .cloned().unwrap_or_else(RecurringTask::new);
        tasks.report_progress    = task_map.get("report_progress")   .cloned().unwrap_or_else(RecurringTask::new);
        tasks
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
            run_initiated:        Instant::now(),
            last_report_progress: Instant::now(),
            tasks_last_performed: Instant::now(),
            call_count_this_run:  0,
        }
    }

    fn report_progress(&mut self, time: f64)
    {
        if self.call_count_this_run > 0 {
            let hours = self.last_report_progress.elapsed().as_secs_f64() / 3600.0;
            println!("");
            println!("\truntime so far ....... {:0.3} hours", self.run_initiated.elapsed().as_secs_f64() / 3600.0);
            println!("");
        }
        self.last_report_progress = Instant::now();
        self.report_progress.advance(10.0);
    }
}
