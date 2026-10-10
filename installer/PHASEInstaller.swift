import AppKit
import Foundation
import SwiftUI

private enum Theme {
    static let background = Color(red: 0.988, green: 0.988, blue: 0.996)
    static let sidebar = Color(red: 0.957, green: 0.965, blue: 0.984)
    static let border = Color(red: 0.894, green: 0.910, blue: 0.941)
    static let ink = Color(red: 0.059, green: 0.078, blue: 0.188)
    static let muted = Color(red: 0.290, green: 0.318, blue: 0.408)
    static let faint = Color(red: 0.549, green: 0.584, blue: 0.722)
    static let blue = Color(red: 0.102, green: 0.310, blue: 0.878)
    static let blueTint = Color(red: 0.910, green: 0.933, blue: 1.000)
    static let green = Color(red: 0.176, green: 0.729, blue: 0.431)
}

@MainActor
final class InstallerModel: ObservableObject {
    let steps = ["Welcome", "MATLAB", "SNAP", "Python", "Destination", "Installation", "Finish"]
    let resources = Bundle.main.resourceURL!
    @Published var page = 0
    @Published var matlab = ""
    @Published var snap = "/Applications/esa-snap/bin/gpt"
    @Published var python = ""
    @Published var destination = FileManager.default.homeDirectoryForCurrentUser
        .appendingPathComponent("Library/Application Support/PHASE").path
    @Published var error = ""
    @Published var details = ""
    @Published var installing = false
    @Published var installed = false
    private var activeProcess: Process?

    var bundledRuntime: Bool {
        FileManager.default.fileExists(atPath: resources.appendingPathComponent("StaMPS/matlab/stamps.m").path)
            && FileManager.default.fileExists(atPath: resources.appendingPathComponent("TRAIN/matlab/aps_linear.m").path)
    }

    init() {
        let apps = (try? FileManager.default.contentsOfDirectory(atPath: "/Applications")) ?? []
        for name in apps.filter({ $0.hasPrefix("MATLAB_R") && $0.hasSuffix(".app") }).sorted().reversed() {
            let candidate = "/Applications/\(name)/bin/matlab"
            if FileManager.default.isExecutableFile(atPath: candidate) {
                matlab = candidate
                break
            }
        }
        let finder = resources.appendingPathComponent("find-python-macos.sh").path
        let process = Process()
        process.executableURL = URL(fileURLWithPath: "/bin/sh")
        process.arguments = [finder]
        let output = Pipe()
        process.standardOutput = output
        process.standardError = Pipe()
        if (try? process.run()) != nil {
            process.waitUntilExit()
            if process.terminationStatus == 0 {
                python = String(data: output.fileHandleForReading.readDataToEndOfFile(),
                                encoding: .utf8)?.trimmingCharacters(in: .whitespacesAndNewlines) ?? ""
            }
        }
    }

    func chooseExecutable(for field: String) {
        let panel = NSOpenPanel()
        panel.canChooseFiles = true
        panel.canChooseDirectories = false
        panel.treatsFilePackagesAsDirectories = true
        panel.allowsOtherFileTypes = true
        panel.message = "Select the \(field) executable"
        guard panel.runModal() == .OK, let path = panel.url?.path else { return }
        switch field {
        case "MATLAB": matlab = path
        case "SNAP gpt": snap = path
        default: python = path
        }
        error = ""
    }

    func chooseDestination() {
        let panel = NSOpenPanel()
        panel.canChooseFiles = false
        panel.canChooseDirectories = true
        panel.canCreateDirectories = true
        panel.message = "Choose the PHASE installation folder"
        if panel.runModal() == .OK, let path = panel.url?.path {
            destination = path
            error = ""
        }
    }

    func next() {
        error = ""
        if page == 1 && !FileManager.default.isExecutableFile(atPath: matlab) {
            error = "Select a valid MATLAB executable to continue."
            return
        }
        if page == 2 && !FileManager.default.isExecutableFile(atPath: snap) {
            error = "Select a valid ESA SNAP gpt executable to continue."
            return
        }
        if page == 3 && !validPython() {
            error = "Select Python 3.10 or newer with venv to continue."
            return
        }
        if page == 4 {
            let normalized = (destination as NSString).expandingTildeInPath
            if normalized.isEmpty || normalized == "/" || normalized == FileManager.default.homeDirectoryForCurrentUser.path {
                error = "Choose a dedicated installation folder."
                return
            }
            destination = normalized
            page = 5
            startInstallation()
            return
        }
        page = min(page + 1, 6)
    }

    private func validPython() -> Bool {
        guard FileManager.default.isExecutableFile(atPath: python) else { return false }
        let process = Process()
        process.executableURL = URL(fileURLWithPath: python)
        process.arguments = ["-c", "import sys,venv; sys.exit(0 if sys.version_info >= (3,10) else 1)"]
        process.standardOutput = Pipe()
        process.standardError = Pipe()
        guard (try? process.run()) != nil else { return false }
        process.waitUntilExit()
        return process.terminationStatus == 0
    }

    func startInstallation() {
        guard !installing else { return }
        installing = true
        installed = false
        error = ""
        details = "Preparing PHASE engine and Python environment…"
        let process = Process()
        process.executableURL = URL(fileURLWithPath: python)
        process.arguments = [resources.appendingPathComponent("install-phase-unix.py").path,
                             "--source", resources.appendingPathComponent("engine").path,
                             "--python", python, "--matlab", matlab, "--gpt", snap,
                             "--prefix", destination]
        if bundledRuntime {
            process.arguments! += ["--stamps", resources.appendingPathComponent("StaMPS").path,
                                   "--train", resources.appendingPathComponent("TRAIN").path]
        }
        let log = FileManager.default.temporaryDirectory
            .appendingPathComponent("phase-install-\(UUID().uuidString).log")
        FileManager.default.createFile(atPath: log.path, contents: nil)
        guard let handle = FileHandle(forWritingAtPath: log.path) else {
            installing = false
            error = "Could not create the installation log."
            return
        }
        process.standardOutput = handle
        process.standardError = handle
        process.terminationHandler = { [weak self] finished in
            try? handle.close()
            let output = (try? String(contentsOf: log, encoding: .utf8)) ?? ""
            DispatchQueue.main.async {
                guard let self else { return }
                self.activeProcess = nil
                self.installing = false
                self.details = output
                if finished.terminationStatus == 0 {
                    self.installed = true
                    self.page = 6
                } else {
                    self.error = output.isEmpty ? "Installation failed. Check the selected paths." : output
                }
            }
        }
        do {
            try process.run()
            activeProcess = process
        } catch {
            try? handle.close()
            installing = false
            self.error = "Unable to start installation: \(error.localizedDescription)"
        }
    }

    func launch() {
        let app = URL(fileURLWithPath: destination).appendingPathComponent("PHASE.app")
        if NSWorkspace.shared.open(app) {
            NSApplication.shared.terminate(nil)
        } else {
            error = "PHASE was installed, but could not launch automatically. Open it from ~/Applications."
        }
    }
}

private struct Card<Content: View>: View {
    @ViewBuilder let content: Content
    var body: some View {
        content
            .padding(22)
            .frame(maxWidth: .infinity, alignment: .leading)
            .background(.white, in: RoundedRectangle(cornerRadius: 12))
            .overlay(RoundedRectangle(cornerRadius: 12).stroke(Theme.border, lineWidth: 1))
    }
}

private struct PhaseButtonStyle: ButtonStyle {
    var primary = false

    func makeBody(configuration: Configuration) -> some View {
        configuration.label
            .font(.system(size: 12, weight: .semibold))
            .foregroundStyle(primary ? Color.white : Theme.blue)
            .padding(.horizontal, 19)
            .frame(minWidth: 100, minHeight: 36)
            .background(primary ? Theme.blue : Color.white,
                        in: RoundedRectangle(cornerRadius: 9))
            .overlay(RoundedRectangle(cornerRadius: 9)
                .stroke(primary ? Theme.blue : Theme.border, lineWidth: 1))
            .opacity(configuration.isPressed ? 0.78 : 1)
            .contentShape(RoundedRectangle(cornerRadius: 9))
    }
}

private struct PathField: View {
    let title: String
    @Binding var value: String
    let browse: () -> Void
    var body: some View {
        VStack(alignment: .leading, spacing: 10) {
            Text(title).font(.system(size: 12, weight: .semibold)).foregroundStyle(Theme.muted)
            HStack(spacing: 10) {
                TextField("Select a path", text: $value)
                    .textFieldStyle(.plain)
                    .font(.system(size: 12, design: .monospaced))
                    .foregroundStyle(Theme.ink)
                    .padding(12)
                    .background(.white, in: RoundedRectangle(cornerRadius: 8))
                    .overlay(RoundedRectangle(cornerRadius: 8).stroke(Theme.border))
                Button("Browse…", action: browse).buttonStyle(PhaseButtonStyle())
            }
        }
    }
}

private struct InstallerView: View {
    @StateObject private var model = InstallerModel()
    private let titles = ["PHASE macOS Installer", "MATLAB", "ESA SNAP", "Python", "Destination folder", "Installation", "Installation complete"]
    private let subtitles = [
        "A single PHASE installation for all your projects.",
        "Choose the MATLAB installation used by PHASE.",
        "Choose the ESA SNAP Graph Processing Tool.",
        "PHASE uses Python for downloads and data utilities.",
        "Choose where the app and processing engine will live.",
        "Installing PHASE and configuring its runtime.",
        "PHASE is ready to launch."
    ]

    var body: some View {
        HStack(spacing: 0) {
            sidebar.frame(width: 278)
            VStack(spacing: 0) {
                header
                ScrollView {
                    VStack(alignment: .leading, spacing: 24) {
                        Text(titles[model.page]).font(.system(size: 29, weight: .light)).foregroundStyle(Theme.ink)
                        Text(subtitles[model.page]).font(.system(size: 14)).foregroundStyle(Theme.muted)
                        pageContent
                        if !model.error.isEmpty {
                            Text(model.error)
                                .font(.system(size: 12))
                                .foregroundStyle(.red)
                                .textSelection(.enabled)
                                .frame(maxWidth: .infinity, alignment: .leading)
                        }
                    }
                    .frame(maxWidth: .infinity, alignment: .leading)
                    .padding(38)
                }
                footer
            }
            .background(Theme.background)
        }
        .frame(width: 980, height: 700)
        .background(Theme.background)
    }

    private var sidebar: some View {
        VStack(alignment: .leading, spacing: 0) {
            if let logo = NSImage(contentsOfFile: model.resources.appendingPathComponent("PHASE_logo.png").path) {
                Image(nsImage: logo).resizable().scaledToFit().frame(height: 62)
                    .frame(maxWidth: .infinity).padding(16)
                    .background(.white, in: RoundedRectangle(cornerRadius: 12))
                    .overlay(RoundedRectangle(cornerRadius: 12).stroke(Theme.border))
                    .padding(.bottom, 32)
            }
            Text("SETUP  ·  PIPELINE")
                .font(.system(size: 10, weight: .semibold, design: .monospaced))
                .foregroundStyle(Theme.faint).padding(.bottom, 16)
            ForEach(model.steps.indices, id: \.self) { index in
                HStack(spacing: 15) {
                    ZStack {
                        Circle().fill(index <= model.page ? Theme.blue : .white)
                        Circle().stroke(index <= model.page ? Theme.blue : Theme.border, lineWidth: 1)
                        Text(index < model.page ? "✓" : String(format: "%02d", index + 1))
                            .font(.system(size: 10, weight: .semibold, design: .monospaced))
                            .foregroundStyle(index <= model.page ? .white : Theme.faint)
                    }.frame(width: 28, height: 28)
                    Text(model.steps[index])
                        .font(.system(size: 13, weight: index == model.page ? .semibold : .regular))
                        .foregroundStyle(index == model.page ? Theme.ink : Theme.muted)
                    Spacer()
                }
                .padding(.vertical, 7)
            }
            Spacer()
            Card {
                VStack(alignment: .leading, spacing: 5) {
                    Text("PERSISTENT SCATTERER").font(.system(size: 10, weight: .bold)).foregroundStyle(Theme.ink)
                    Text("Highly Automated Suite for Environmental Monitoring")
                        .font(.system(size: 10)).foregroundStyle(Theme.faint)
                }
            }.padding(.bottom, 14)
            Text("PHASE  v7.0.0 preview  ·  Roberto Monti · pyccino")
                .font(.system(size: 9, design: .monospaced)).foregroundStyle(Theme.faint)
        }
        .padding(24)
        .background(Theme.sidebar)
        .overlay(alignment: .trailing) { Theme.border.frame(width: 1) }
    }

    private var header: some View {
        HStack {
            Text("STEP \(model.page + 1)  ·  \(model.steps[model.page].uppercased())")
                .foregroundStyle(Theme.blue)
            Spacer()
            Text("// phase-installer // Roberto Monti · pyccino")
                .foregroundStyle(Theme.faint)
        }
        .font(.system(size: 10, weight: .semibold, design: .monospaced))
        .padding(.horizontal, 38).frame(height: 52)
        .background(.white)
        .overlay(alignment: .bottom) { Theme.border.frame(height: 1) }
    }

    @ViewBuilder private var pageContent: some View {
        switch model.page {
        case 0:
            Card {
                VStack(alignment: .leading, spacing: 14) {
                    Text("What this wizard does").font(.system(size: 15, weight: .semibold))
                    Label("Verify MATLAB, SNAP and Python", systemImage: "checkmark.circle")
                    Label("Install the PHASE app and processing engine", systemImage: "square.and.arrow.down")
                    Label(model.bundledRuntime ? "Include Apple Silicon StaMPS and TRAIN" : "Install without PSI runtime", systemImage: "waveform.path")
                    Label("Create PHASE.app in your user Applications folder", systemImage: "app")
                }.font(.system(size: 13)).foregroundStyle(Theme.ink)
            }
            Card {
                VStack(alignment: .leading, spacing: 6) {
                    Text("Install once. Create projects anywhere.").font(.system(size: 14, weight: .semibold))
                    Text("Your project folders and outputs stay separate from the app, including after updates.")
                        .font(.system(size: 12))
                }.foregroundStyle(Theme.ink)
            }.background(Theme.blueTint)
        case 1:
            Card {
                VStack(alignment: .leading, spacing: 16) {
                    PathField(title: "MATLAB executable", value: $model.matlab) { model.chooseExecutable(for: "MATLAB") }
                    Text(model.matlab.isEmpty ? "MATLAB was not detected. Select it manually." : "Detected MATLAB. You can choose another installation.")
                        .font(.system(size: 12)).foregroundStyle(Theme.muted)
                }
            }
        case 2:
            Card {
                VStack(alignment: .leading, spacing: 16) {
                    PathField(title: "SNAP gpt executable", value: $model.snap) { model.chooseExecutable(for: "SNAP gpt") }
                    Text("ESA SNAP must be installed separately; PHASE will use this path for preprocessing.")
                        .font(.system(size: 12)).foregroundStyle(Theme.muted)
                }
            }
        case 3:
            Card {
                VStack(alignment: .leading, spacing: 16) {
                    PathField(title: "Python 3.10+ executable", value: $model.python) { model.chooseExecutable(for: "Python") }
                    Text("The installer creates an isolated Python environment for PHASE.")
                        .font(.system(size: 12)).foregroundStyle(Theme.muted)
                }
            }
        case 4:
            Card {
                VStack(alignment: .leading, spacing: 17) {
                    PathField(title: "PHASE installation folder", value: $model.destination) { model.chooseDestination() }
                    Text("One managed installation. PHASE.app will also appear in ~/Applications. Projects may be stored anywhere and are not changed by reinstalling or updating PHASE.")
                        .font(.system(size: 12)).foregroundStyle(Theme.muted)
                    Text(model.bundledRuntime ? "StaMPS PSI runtime included" : "StaMPS PSI runtime not included in this installer")
                        .font(.system(size: 12, weight: .medium))
                        .foregroundStyle(model.bundledRuntime ? Theme.green : .orange)
                }
            }
        case 5:
            Card {
                VStack(alignment: .leading, spacing: 20) {
                    if model.installing { ProgressView().controlSize(.large) }
                    Text(model.installing ? "Installing PHASE…" : "Installation stopped")
                        .font(.system(size: 17, weight: .semibold))
                    Text(model.installing ? "Copying the app and native runtime, then preparing Python packages. This may take a few minutes." : "Check the error below and retry.")
                        .font(.system(size: 12)).foregroundStyle(Theme.muted)
                    if !model.installing {
                        Button("Retry installation") { model.startInstallation() }
                            .buttonStyle(PhaseButtonStyle(primary: true))
                    }
                }
            }
        default:
            Card {
                VStack(alignment: .leading, spacing: 15) {
                    Label("PHASE is ready", systemImage: "checkmark.circle.fill")
                        .font(.system(size: 19, weight: .medium)).foregroundStyle(Theme.green)
                    Text("The app is in ~/Applications. Its managed files are at:")
                        .font(.system(size: 12)).foregroundStyle(Theme.muted)
                    Text(model.destination).font(.system(size: 12, design: .monospaced)).textSelection(.enabled)
                    Text("Create or open a project from PHASE; project folders can be anywhere.")
                        .font(.system(size: 12)).foregroundStyle(Theme.muted)
                }
            }
        }
    }

    private var footer: some View {
        HStack {
            Button("Back") { model.page -= 1 }
                .buttonStyle(PhaseButtonStyle())
                .disabled(model.page == 0 || model.page >= 5)
            Spacer()
            if model.page == 6 {
                Button("Close") { NSApplication.shared.terminate(nil) }
                    .buttonStyle(PhaseButtonStyle())
                Button("Launch PHASE") { model.launch() }
                    .buttonStyle(PhaseButtonStyle(primary: true))
            } else if model.page < 5 {
                Button(model.page == 4 ? "Install PHASE" : "Next") { model.next() }
                    .buttonStyle(PhaseButtonStyle(primary: true))
            }
        }
        .padding(.horizontal, 38).frame(height: 76)
        .background(.white)
        .overlay(alignment: .top) { Theme.border.frame(height: 1) }
    }
}

@main
struct PHASEInstallerApp: App {
    var body: some Scene {
        WindowGroup {
            InstallerView()
        }
        .windowStyle(.titleBar)
        .windowResizability(.contentSize)
    }
}
