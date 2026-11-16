//! Circular delay history buffer with cubic Hermite spline interpolation for Delay-Differential Equations (DDEs).

/// An entry in the temporal history buffer.
#[derive(Debug, Clone, Copy)]
pub struct HistoryPoint {
    /// Time stamp in minutes.
    pub time: f64,
    /// State value at this time.
    pub value: f64,
    /// Time derivative dy/dt at this time.
    pub derivative: f64,
}

/// Ring-buffer storing continuous historical trajectories for delayed physiological processes.
#[derive(Debug, Clone)]
pub struct DelayBuffer {
    buffer: Vec<HistoryPoint>,
    capacity: usize,
    head: usize,
    count: usize,
}

impl DelayBuffer {
    /// Create a new delay buffer with a fixed storage capacity.
    pub fn new(capacity: usize) -> Self {
        Self {
            buffer: vec![
                HistoryPoint {
                    time: 0.0,
                    value: 0.0,
                    derivative: 0.0
                };
                capacity.max(8)
            ],
            capacity: capacity.max(8),
            head: 0,
            count: 0,
        }
    }

    /// Push a new observed time point into the history buffer.
    pub fn push(&mut self, time: f64, value: f64, derivative: f64) {
        let pt = HistoryPoint {
            time,
            value,
            derivative,
        };
        self.buffer[self.head] = pt;
        self.head = (self.head + 1) % self.capacity;
        if self.count < self.capacity {
            self.count += 1;
        }
    }

    /// Query the historical value at $t - \tau$ using cubic Hermite spline interpolation.
    ///
    /// If requested time is prior to the oldest recorded sample, clamps to the oldest known value.
    /// If requested time is after the newest sample, extrapolates linearly with the latest derivative.
    pub fn interpolate_delayed(&self, query_time: f64) -> f64 {
        if self.count == 0 {
            return 0.0;
        }

        let newest_idx = (self.head + self.capacity - 1) % self.capacity;
        let newest = self.buffer[newest_idx];

        // If query is after newest point, extrapolate forward
        if query_time >= newest.time {
            let dt = query_time - newest.time;
            return newest.value + dt * newest.derivative;
        }

        let oldest_idx = if self.count < self.capacity {
            0
        } else {
            self.head
        };
        let oldest = self.buffer[oldest_idx];

        // If query is older than buffer history, clamp to oldest
        if query_time <= oldest.time {
            return oldest.value;
        }

        // Binary search or linear scan through circular buffer
        // Since count is moderate (e.g. 500-2000 points), scan from newest backward
        let mut right_pt = newest;
        for i in 1..self.count {
            let idx = (self.head + self.capacity - 1 - i) % self.capacity;
            let left_pt = self.buffer[idx];

            if left_pt.time <= query_time && query_time <= right_pt.time {
                let h = right_pt.time - left_pt.time;
                if h.abs() < 1e-12 {
                    return left_pt.value;
                }

                // Normalized local coordinate u in [0, 1]
                let u = (query_time - left_pt.time) / h;
                let u2 = u * u;
                let u3 = u2 * u;

                // Hermite basis polynomials
                let h00 = 2.0 * u3 - 3.0 * u2 + 1.0;
                let h10 = u3 - 2.0 * u2 + u;
                let h01 = -2.0 * u3 + 3.0 * u2;
                let h11 = u3 - u2;

                return h00 * left_pt.value
                    + h10 * h * left_pt.derivative
                    + h01 * right_pt.value
                    + h11 * h * right_pt.derivative;
            }

            right_pt = left_pt;
        }

        oldest.value
    }

    /// Total number of stored temporal samples.
    pub fn len(&self) -> usize {
        self.count
    }

    /// Returns true if no samples have been recorded.
    pub fn is_empty(&self) -> bool {
        self.count == 0
    }

    /// Clear all recorded history.
    pub fn clear(&mut self) {
        self.head = 0;
        self.count = 0;
    }
}
