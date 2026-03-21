import { Routes, Route } from 'react-router-dom';
import { Dashboard } from './Dashboard';

export default function App() {
  return (
    <Routes>
      <Route path="/:projectName/:docId" element={<Dashboard />} />
      <Route path="/:projectName" element={<Dashboard />} />
      <Route path="/" element={<Dashboard />} />
    </Routes>
  );
}
